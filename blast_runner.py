"""
PRNP-OrthoMiner — Motor del pipeline BLAST
Adaptación a Python del script BASH original de Ana Rosa Cortazar (marzo 2021).

Flujo por cada archivo .gz:
  1. Descomprimir el ensamblaje genómico
  2. Ejecutar blastn (genoma como query, PrP como subject)
  3. Localizar el contig con alineamiento significativo
  4. Extraer ese contig y generar _output.txt + _outputFA.fasta
  5. Limpiar archivos temporales
"""

import os
import gzip
import shutil
import subprocess
import glob
import queue
import threading
import uuid

# -------------------------------------------------------------------
# Registro de jobs en curso  {job_id: queue.Queue}
# -------------------------------------------------------------------
_jobs: dict = {}


# -------------------------------------------------------------------
# Utilidades de ficheros
# -------------------------------------------------------------------

def find_gz_files(base_path: str, species: str) -> list:
    """Busca archivos .gz en <base_path>/<species>/ y sus subdirectorios."""
    species_dir = os.path.join(base_path.rstrip('/\\'), species)
    if not os.path.isdir(species_dir):
        return []
    files = glob.glob(os.path.join(species_dir, '*.gz'))
    for entry in os.scandir(species_dir):
        if entry.is_dir():
            files.extend(glob.glob(os.path.join(entry.path, '*.gz')))
    return sorted(set(files))


def decompress_gz(gz_path: str) -> str:
    """Descomprime <fichero>.gz -> <fichero>. Devuelve la ruta descomprimida."""
    out_path = gz_path[:-3]  # quita .gz
    with gzip.open(gz_path, 'rb') as f_in:
        with open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return out_path


def write_reference_fasta(name: str, sequence: str, path: str) -> None:
    """Escribe la secuencia PrP de referencia en formato FASTA."""
    with open(path, 'w') as f:
        f.write(f'>{name}\n{sequence}\n')


def _cleanup(*paths):
    """Elimina archivos temporales ignorando errores."""
    for p in paths:
        try:
            if p and os.path.exists(p):
                os.remove(p)
        except OSError:
            pass


# -------------------------------------------------------------------
# BLAST
# -------------------------------------------------------------------

def run_blastn(genome_path: str, ref_fasta: str, out_path: str) -> tuple:
    """
    Ejecuta: blastn -query <genoma> -subject <PrP.fa> -outfmt 1 -out <salida>
    Devuelve (returncode, stderr).
    """
    cmd = [
        'blastn',
        '-query',   genome_path,
        '-subject', ref_fasta,
        '-outfmt',  '1',
        '-out',     out_path
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        return result.returncode, result.stderr
    except FileNotFoundError:
        raise RuntimeError(
            "BLAST+ no encontrado. Instala BLAST+ y asegúrate de que "
            "'blastn' está en el PATH del sistema."
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError("Tiempo límite agotado ejecutando BLAST (máx. 2 horas).")


# -------------------------------------------------------------------
# Parseo de resultados y extracción de secuencia
# -------------------------------------------------------------------

def parse_and_extract(
    blast_out: str,
    genome_path: str,
    output_dir: str,
    base_name: str
) -> bool:
    """
    Analiza la salida BLAST (-outfmt 1).
    - Si hay alineamiento significativo genera _output.txt y _outputFA.fasta.
    - Devuelve True si se encontró alineamiento, False si no.
    """
    with open(blast_out, 'r') as f:
        lines = f.readlines()

    # Buscar línea "Sequences producing significant alignments:"
    target_idx = None
    for i, line in enumerate(lines):
        if 'Sequences producing significant alignments:' in line:
            target_idx = i
            break

    if target_idx is None:
        return False

    # 5 líneas antes está la cabecera "Query= <contig>"
    start_idx = max(0, target_idx - 5)

    # --- Generar _output.txt ---
    out_txt = os.path.join(output_dir, base_name + '_output.txt')
    with open(out_txt, 'w') as f:
        # Cabecera (incluye la línea Query=)
        for line in lines[start_idx:target_idx]:
            f.write(line)
        # Sección de alineamientos hasta el próximo "Query="
        writing = False
        for i, line in enumerate(lines):
            if 'Sequences producing significant alignments:' in line:
                writing = True
            if writing:
                if line.startswith('Query=') and i > target_idx:
                    break
                f.write(line)

    # Extraer nombre del contig desde la línea Query=
    buscar = None
    header_line = lines[start_idx] if start_idx < len(lines) else ''
    if 'Query= ' in header_line:
        buscar = header_line.split('Query= ', 1)[1].strip()

    # --- Generar _outputFA.fasta ---
    out_fasta = os.path.join(output_dir, base_name + '_outputFA.fasta')
    with open(out_fasta, 'w') as ff:
        ff.write('>Secuencia_Extendida\n')
        if buscar:
            with open(genome_path, 'r') as gf:
                in_target = False
                for line in gf:
                    if line.startswith('>'):
                        if in_target:
                            break
                        if buscar in line:
                            in_target = True
                            ff.write(line)
                    elif in_target:
                        ff.write(line)

    return True


# -------------------------------------------------------------------
# Pipeline completo (generador de eventos de progreso)
# -------------------------------------------------------------------

def _pipeline_generator(base_path, species, ref_name, ref_sequence, output_path):
    """
    Generador que ejecuta el pipeline paso a paso y emite tuplas
    (event_type, message) donde event_type es:
      'info', 'progress', 'step', 'success', 'warning', 'error', 'done'
    """
    yield 'info',     'Iniciando PRNP-OrthoMiner pipeline'
    yield 'info',     f'Especie: {species}'
    yield 'info',     f'Referencia PrP: {ref_name}'
    yield 'info',     f'Ruta base: {base_path}'
    yield 'info',     f'Ruta salida: {output_path}'

    if not os.path.isdir(base_path):
        yield 'error', f'Ruta base no encontrada: {base_path}'
        return

    gz_files = find_gz_files(base_path, species)
    if not gz_files:
        yield 'error', f'Sin archivos .gz en: {os.path.join(base_path, species)}'
        return

    yield 'info', f'Archivos .gz encontrados: {len(gz_files)}'
    os.makedirs(output_path, exist_ok=True)

    # Escribir FASTA de referencia
    ref_fasta = os.path.join(output_path, '_prp_ref_temp.fa')
    write_reference_fasta(ref_name, ref_sequence, ref_fasta)

    found = 0
    errors = 0

    for i, gz_file in enumerate(gz_files):
        filename  = os.path.basename(gz_file)
        base_name = filename[:-3] if filename.endswith('.gz') else filename
        yield 'progress', f'[{i+1}/{len(gz_files)}] {filename}'

        # 1. Descomprimir
        try:
            yield 'step', 'Descomprimiendo…'
            genome_path = decompress_gz(gz_file)
        except Exception as e:
            yield 'error', f'Error descomprimiendo {filename}: {e}'
            errors += 1
            continue

        # 2. BLAST
        blast_out = genome_path + '_outputTemp'
        try:
            yield 'step', 'Ejecutando BLAST…'
            rc, stderr = run_blastn(genome_path, ref_fasta, blast_out)
            if rc != 0:
                yield 'error', f'BLAST falló en {filename}: {stderr[:300]}'
                errors += 1
                _cleanup(genome_path, blast_out)
                continue
        except RuntimeError as e:
            yield 'error', str(e)
            _cleanup(genome_path, blast_out, ref_fasta)
            return

        # 3. Parseo y extracción
        try:
            yield 'step', 'Procesando resultados…'
            hit = parse_and_extract(blast_out, genome_path, output_path, base_name)
            if hit:
                found += 1
                yield 'success', (
                    f'¡Alineamiento encontrado! → '
                    f'{base_name}_output.txt  |  {base_name}_outputFA.fasta'
                )
            else:
                yield 'warning', 'Sin alineamientos significativos'
        except Exception as e:
            yield 'error', f'Error procesando {filename}: {e}'
            errors += 1

        _cleanup(genome_path, blast_out)

    _cleanup(ref_fasta)
    yield 'done', (
        f'Pipeline completado — '
        f'Encontrados: {found}/{len(gz_files)} | Errores: {errors}'
    )


# -------------------------------------------------------------------
# API pública: jobs asíncronos con cola
# -------------------------------------------------------------------

def start_job(base_path, species, ref_name, ref_sequence, output_path) -> str:
    """Lanza el pipeline en un hilo background. Devuelve el job_id."""
    job_id = uuid.uuid4().hex[:8]
    q: queue.Queue = queue.Queue()
    _jobs[job_id] = q

    def worker():
        try:
            for etype, msg in _pipeline_generator(
                base_path, species, ref_name, ref_sequence, output_path
            ):
                q.put({'type': etype, 'message': msg})
        except Exception as e:
            q.put({'type': 'error', 'message': str(e)})
        finally:
            q.put(None)  # sentinel — fin del stream

    threading.Thread(target=worker, daemon=True).start()
    return job_id


def get_job_queue(job_id: str):
    """Devuelve la cola del job o None si no existe."""
    return _jobs.get(job_id)


def cleanup_job(job_id: str):
    """Elimina el job de la tabla."""
    _jobs.pop(job_id, None)
