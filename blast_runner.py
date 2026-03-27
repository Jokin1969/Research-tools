"""
PRNP-OrthoMiner — Motor del pipeline BLAST
Adaptación mejorada del script BASH original de Ana Rosa Cortazar (marzo 2021).

Mejoras respecto al original:
  - Extrae coordenadas exactas del alineamiento (qstart/qend) via -outfmt 6
  - La extensión ±200 bp se aplica sobre coordenadas reales, no con sed fragile
  - Selecciona el mejor hit por bitscore si hay varios contigs con alineamiento
  - El informe de texto (_output.txt) es fiel al formato original (-outfmt 1)
  - El FASTA de salida (_outputFA.fasta) incluye la región extendida con cabecera informativa

Flujo por cada archivo .gz:
  1. Descomprimir el ensamblaje genómico
  2. blastn -outfmt 1  -> informe legible (_output.txt)
  3. blastn -outfmt 6  -> coordenadas precisas del mejor alineamiento
  4. Extraer contig completo del FASTA y recortar qstart-200 .. qend+200
  5. Guardar _outputFA.fasta y limpiar temporales
"""

import os
import gzip
import shutil
import subprocess
import glob
import queue
import threading
import uuid
from pathlib import Path

FLANK_BP = 200          # Nucleotídos de extensión en cada extremo

# ---------------------------------------------------------------------------
# Registro de jobs  {job_id: queue.Queue}
# ---------------------------------------------------------------------------
_jobs: dict = {}

# Archivos temporales conservados cuando no hay hit (para reintento)
# {job_id: {'work_dir': str, 'species': str}}
_kept_files: dict = {}


def get_kept_files(job_id: str):
    return _kept_files.get(job_id)


def consume_kept_files(job_id: str):
    """Quita y devuelve la entrada sin borrar los archivos (para reintento)."""
    return _kept_files.pop(job_id, None)


def delete_kept_files(job_id: str) -> bool:
    entry = _kept_files.pop(job_id, None)
    if entry and os.path.isdir(entry['work_dir']):
        shutil.rmtree(entry['work_dir'], ignore_errors=True)
        return True
    return False


# ---------------------------------------------------------------------------
# Utilidades de ficheros
# ---------------------------------------------------------------------------

def find_gz_files(base_path: str, species: str) -> list:
    """
    Busca archivos .gz en:  <base_path>/<species>/*.gz
    y en subdirectorios directos de esa carpeta (como el original: WGS2/*/*.gz).
    """
    species_dir = os.path.join(base_path.rstrip('/\\'), species)
    if not os.path.isdir(species_dir):
        return []
    files = glob.glob(os.path.join(species_dir, '*.gz'))
    for entry in os.scandir(species_dir):
        if entry.is_dir():
            files.extend(glob.glob(os.path.join(entry.path, '*.gz')))
    return sorted(set(files))


def decompress_gz(gz_path: str) -> str:
    """Descomprime <archivo>.gz y devuelve la ruta del archivo descomprimido."""
    out_path = gz_path[:-3]
    with gzip.open(gz_path, 'rb') as f_in:
        with open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    return out_path


def write_reference_fasta(name: str, sequence: str, path: str) -> None:
    """Escribe la secuencia PrP de referencia en formato FASTA."""
    with open(path, 'w') as f:
        f.write(f'>{name}\n{sequence}\n')


def _cleanup(*paths):
    for p in paths:
        try:
            if p and os.path.exists(p):
                os.remove(p)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# BLAST
# ---------------------------------------------------------------------------

def _blastn(genome: str, ref: str, outfmt: str, out_path: str) -> tuple:
    """
    Ejecuta blastn con el formato indicado.
    Devuelve (returncode, stderr).
    Lanza RuntimeError si BLAST+ no está instalado o hay timeout.
    """
    cmd = ['blastn', '-query', genome, '-subject', ref,
           '-outfmt', outfmt, '-out', out_path]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        return r.returncode, r.stderr
    except FileNotFoundError:
        raise RuntimeError(
            "BLAST+ no encontrado. Instala NCBI BLAST+ y asegúrate de que "
            "'blastn' está en el PATH del sistema."
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            "Tiempo límite agotado ejecutando BLAST (máx. 2 horas)."
        )


def run_blast_dual(genome: str, ref: str, txt_out: str, tab_out: str) -> tuple:
    """
    Ejecuta BLAST dos veces:
      1) -outfmt 1  -> informe de texto legible (compatible con el original)
      2) -outfmt 6  -> tabla de coordenadas para extracción precisa
    Devuelve (returncode, stderr) del último paso.
    """
    rc, err = _blastn(genome, ref, '1', txt_out)
    if rc != 0:
        return rc, err
    rc, err = _blastn(genome, ref, '6', tab_out)
    return rc, err


# ---------------------------------------------------------------------------
# Parseo de salida tabular (-outfmt 6)
# ---------------------------------------------------------------------------
# Columnas estándar de -outfmt 6:
#   0 qseqid  1 sseqid  2 pident  3 length  4 mismatch  5 gapopen
#   6 qstart  7 qend    8 sstart  9 send   10 evalue   11 bitscore

def best_hit_from_tabular(tab_path: str):
    """
    Lee el fichero tabular (-outfmt 6) y devuelve el mejor alineamiento
    (el de mayor bitscore), o None si no hay hits.
    """
    best = None
    if not os.path.exists(tab_path):
        return None
    with open(tab_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 12:
                continue
            hit = {
                'qseqid':   parts[0],
                'sseqid':   parts[1],
                'pident':   float(parts[2]),
                'length':   int(parts[3]),
                'qstart':   int(parts[6]),   # 1-based
                'qend':     int(parts[7]),   # 1-based
                'sstart':   int(parts[8]),
                'send':     int(parts[9]),
                'evalue':   float(parts[10]),
                'bitscore': float(parts[11])
            }
            hit['strand'] = '+' if hit['sstart'] <= hit['send'] else '-'
            if best is None or hit['bitscore'] > best['bitscore']:
                best = hit
    return best


# ---------------------------------------------------------------------------
# Parseo del informe de texto (-outfmt 1) para _output.txt
# (fiel al comportamiento original del script BASH)
# ---------------------------------------------------------------------------

def extract_text_report(txt_blast: str, out_txt: str) -> bool:
    """
    Replica la lógica del BASH original:
      - Localiza la línea 'Sequences producing significant alignments:'
      - Extrae 5 líneas de cabecera (Query= incluido)
      - Extrae la sección de alineamientos hasta el siguiente 'Query='
      - Escribe el resultado en out_txt
    Devuelve True si se encontró alineamiento, False si no.
    """
    with open(txt_blast) as f:
        lines = f.readlines()

    target_idx = None
    for i, line in enumerate(lines):
        if 'Sequences producing significant alignments:' in line:
            target_idx = i
            break

    if target_idx is None:
        return False

    start_idx = max(0, target_idx - 5)

    with open(out_txt, 'w') as out:
        # Cabecera (incluye línea Query=)
        for line in lines[start_idx:target_idx]:
            out.write(line)
        # Sección de alineamientos hasta el siguiente bloque Query=
        writing = False
        for i, line in enumerate(lines):
            if 'Sequences producing significant alignments:' in line:
                writing = True
            if writing:
                if line.startswith('Query=') and i > target_idx:
                    break
                out.write(line)

    return True


# ---------------------------------------------------------------------------
# Extracción de secuencia FASTA con extensión de flanqueantes
# ---------------------------------------------------------------------------

def read_fasta_sequences(fasta_path: str) -> dict:
    """
    Lee un FASTA multi-contig y devuelve {seq_id: sequence_string}.
    El seq_id es la primera palabra tras '>' (sin espacios ni \n).
    """
    seqs = {}
    current_id = None
    buffer = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_id is not None:
                    seqs[current_id] = ''.join(buffer)
                current_id = line[1:].split()[0]   # primera palabra del header
                buffer = []
            else:
                buffer.append(line)
    if current_id is not None:
        seqs[current_id] = ''.join(buffer)
    return seqs


def extract_extended_fasta(genome_path: str, hit: dict, out_fasta: str, species: str = ''):
    """
    Extrae el locus candidato de PRNP con ±FLANK_BP bp de contexto flanqueante.

    Usa las coordenadas qstart/qend del mejor hit tabular (base 1).
    Si el contig es más corto que la ventana solicitada se ajusta al extremo.

    Formato de salida (mixed-case):
      - Flancos (±FLANK_BP bp): minúsculas
      - Región alineada (ORF candidato): MAYÚSCULAS

    Devuelve un dict con 'fasta_content', 'genomic_info' y metadatos,
    o None si no se pudo localizar el contig en el genoma.
    """
    seqs = read_fasta_sequences(genome_path)

    # Buscar el contig que coincide con qseqid (coincidencia exacta o parcial)
    contig_id  = hit['qseqid']
    target_seq = seqs.get(contig_id)
    if target_seq is None:
        # Búsqueda parcial por si el ID fue truncado por BLAST
        for sid, seq in seqs.items():
            if contig_id in sid or sid in contig_id:
                target_seq = seq
                contig_id  = sid
                break

    if target_seq is None:
        return None

    # Convertir coordenadas BLAST (base 1) a índices Python (base 0)
    qstart = hit['qstart'] - 1   # inclusive
    qend   = hit['qend']         # exclusive
    if qstart > qend:            # alineamiento en strand negativo: invertir
        qstart, qend = qend - 1, qstart + 1

    # Aplicar extensión flanqueante
    ext_start = max(0, qstart - FLANK_BP)
    ext_end   = min(len(target_seq), qend + FLANK_BP)

    # Mixed-case: flancos en minúsculas, región alineada en MAYÚSCULAS
    upstream   = target_seq[ext_start:qstart].lower()
    orf_region = target_seq[qstart:qend].upper()
    downstream = target_seq[qend:ext_end].lower()
    region     = upstream + orf_region + downstream

    # Cabecera FASTA con nombre de especie
    header_name = f"{species}_PRNP_Extended_Region" if species else "Locus_PRNP_candidato"
    header = (
        f">{header_name} | "
        f"contig={contig_id} | "
        f"coords={ext_start+1}-{ext_end} | "
        f"aln_coords={hit['qstart']}-{hit['qend']} | "
        f"strand={hit.get('strand', '?')} | "
        f"pident={hit['pident']:.1f}% | "
        f"evalue={hit['evalue']} | "
        f"ext=\u00b1{FLANK_BP}bp"
    )

    # Formatear la secuencia en líneas de 80 caracteres (estándar FASTA)
    fasta_lines = [region[i:i+80] for i in range(0, len(region), 80)]
    fasta_content = header + '\n' + '\n'.join(fasta_lines) + '\n'

    # Información genómica para mostrar en la interfaz
    genomic_info = (
        f"Especie        : {species or 'Desconocida'}\n"
        f"Scaffold/Contig: {contig_id}\n"
        f"Regi\u00f3n extendida: {ext_start+1}\u2013{ext_end} (base 1)\n"
        f"Regi\u00f3n ORF     : {hit['qstart']}\u2013{hit['qend']} (coords BLAST, base 1)\n"
        f"Hebra          : {hit.get('strand', '?')}\n"
        f"Identidad      : {hit['pident']:.1f}%\n"
        f"E-value        : {hit['evalue']}\n"
        f"Longitud aln.  : {hit['length']} bp\n"
        f"Extensi\u00f3n      : \u00b1{FLANK_BP} bp\n"
        f"Upstream (min.): {len(upstream)} bp\n"
        f"ORF (MAY\u00daS.)  : {len(orf_region)} bp\n"
        f"Downstream(min): {len(downstream)} bp"
    )

    with open(out_fasta, 'w') as f:
        f.write(fasta_content)

    return {
        'contig':        contig_id,
        'ext_start':     ext_start + 1,
        'ext_end':       ext_end,
        'aln_start':     hit['qstart'],
        'aln_end':       hit['qend'],
        'strand':        hit.get('strand', '?'),
        'pident':        hit['pident'],
        'evalue':        hit['evalue'],
        'fasta_content': fasta_content,
        'genomic_info':  genomic_info
    }


# ---------------------------------------------------------------------------
# Pipeline completo (generador de eventos de progreso)
# ---------------------------------------------------------------------------

def _pipeline_generator(base_path, species, ref_name, ref_sequence, output_path):
    """
    Generador que ejecuta el pipeline paso a paso.
    Emite tuplas (event_type, message) con tipos:
      'info', 'progress', 'step', 'success', 'warning', 'error', 'done'
    """
    yield 'info',  'Iniciando PRNP-OrthoMiner pipeline'
    yield 'info',  f'Especie        : {species}'
    yield 'info',  f'Referencia PrP : {ref_name}'
    yield 'info',  f'Ruta base      : {base_path}'
    yield 'info',  f'Ruta salida    : {output_path}'
    yield 'info',  f'Extensión flanq.: ±{FLANK_BP} bp'

    species_path = os.path.join(base_path.rstrip('/\\'), species)

    if not os.path.isdir(base_path):
        yield 'error', (
            f'Ruta base no encontrada: {base_path} — '
            f'Recuerda que el análisis corre en el servidor; '
            f'las rutas deben ser accesibles desde el servidor, no desde tu máquina local.'
        )
        return

    if not os.path.isdir(species_path):
        yield 'error', f'Directorio de especie no encontrado: {species_path}'
        return

    gz_files = find_gz_files(base_path, species)
    if not gz_files:
        yield 'error', f'Sin archivos .gz en: {species_path}'
        return

    total = len(gz_files)
    yield 'info', f'Archivos .gz encontrados: {total}'
    os.makedirs(output_path, exist_ok=True)

    # Escribir FASTA de referencia
    ref_fasta = os.path.join(output_path, '_prp_ref_temp.fa')
    write_reference_fasta(ref_name, ref_sequence, ref_fasta)

    found = 0
    no_hit = 0
    errors = 0

    for i, gz_file in enumerate(gz_files):
        filename  = os.path.basename(gz_file)
        base_name = filename[:-3] if filename.endswith('.gz') else filename
        yield 'progress', f'[{i+1}/{total}]  {filename}'

        # 1. Descomprimir ---------------------------------------------------
        try:
            yield 'step', 'Descomprimiendo…'
            genome_path = decompress_gz(gz_file)
        except Exception as e:
            yield 'error', f'Error descomprimiendo {filename}: {e}'
            errors += 1
            continue

        # Rutas de salida para este archivo
        txt_blast = genome_path + '_blast_fmt1.tmp'
        tab_blast = genome_path + '_blast_fmt6.tmp'
        out_txt   = os.path.join(output_path, base_name + '_output.txt')
        out_fasta = os.path.join(output_path, base_name + '_outputFA.fasta')

        # 2. BLAST (x2) -----------------------------------------------------
        try:
            yield 'step', 'Ejecutando BLAST (informe de texto + tabla de coordenadas)…'
            rc, stderr = run_blast_dual(genome_path, ref_fasta, txt_blast, tab_blast)
            if rc != 0:
                yield 'error', f'BLAST falló ({filename}): {stderr[:300]}'
                errors += 1
                _cleanup(genome_path, txt_blast, tab_blast)
                continue
        except RuntimeError as e:
            yield 'error', str(e)
            _cleanup(genome_path, txt_blast, tab_blast, ref_fasta)
            return

        # 3. Parseo del informe de texto ------------------------------------
        try:
            yield 'step', 'Analizando informe de texto…'
            has_hit = extract_text_report(txt_blast, out_txt)
        except Exception as e:
            yield 'error', f'Error procesando informe ({filename}): {e}'
            has_hit = False

        if not has_hit:
            yield 'warning', 'Sin alineamientos significativos'
            no_hit += 1
            _cleanup(genome_path, txt_blast, tab_blast)
            continue

        # 4. Mejor hit + extracción de secuencia con flanqueantes -----------
        try:
            yield 'step', f'Extrayendo locus candidato (±{FLANK_BP} bp)…'
            hit = best_hit_from_tabular(tab_blast)
            if hit is None:
                yield 'warning', 'No se leyeron coordenadas del output tabular; FASTA no generado.'
                no_hit += 1
            else:
                result = extract_extended_fasta(genome_path, hit, out_fasta, species)
                if result:
                    found += 1
                    yield 'success', (
                        f'¡Locus PRNP encontrado! '
                        f'contig={hit["qseqid"]} | '
                        f'coords={hit["qstart"]}-{hit["qend"]} | '
                        f'identidad={hit["pident"]:.1f}% | '
                        f'e-value={hit["evalue"]}'
                    )
                    yield 'result', {
                        'species':       species,
                        'genome_file':   base_name,
                        'fasta_content': result['fasta_content'],
                        'genomic_info':  result['genomic_info']
                    }
                    _cleanup(genome_path, txt_blast, tab_blast)
                    yield 'info', 'Análisis detenido: locus encontrado. Archivos temporales eliminados.'
                    break   # parar al primer hit
                else:
                    yield 'warning', 'No se pudo extraer la secuencia del contig del genoma.'
                    no_hit += 1
        except Exception as e:
            yield 'error', f'Error extrayendo secuencia ({filename}): {e}'
            errors += 1

        _cleanup(genome_path, txt_blast, tab_blast)

    # --- Fin del loop -------------------------------------------------------
    _cleanup(ref_fasta)
    yield 'done', (
        f'Pipeline completado — '
        f'Encontrados: {found}/{total} | '
        f'Sin hit: {no_hit}/{total} | '
        f'Errores: {errors}/{total}'
    )


# ---------------------------------------------------------------------------
# API pública: jobs asíncronos con cola
# ---------------------------------------------------------------------------

def start_job(base_path, species, ref_name, ref_sequence, output_path,
              cleanup_dir=None) -> str:
    """Lanza el pipeline en un hilo background y devuelve el job_id.

    Si se pasa cleanup_dir, ese directorio se borra automáticamente
    al finalizar el pipeline (usado en modo subida de archivos).
    """
    job_id = uuid.uuid4().hex[:8]
    q: queue.Queue = queue.Queue()
    _jobs[job_id] = q

    def worker():
        hit_found = False
        try:
            for etype, msg in _pipeline_generator(
                base_path, species, ref_name, ref_sequence, output_path
            ):
                q.put({'type': etype, 'message': msg})
                if etype == 'result':
                    hit_found = True
        except Exception as e:
            q.put({'type': 'error', 'message': str(e)})
        finally:
            q.put(None)   # sentinel — fin del stream
            if cleanup_dir:
                if hit_found:
                    # Hit encontrado: los archivos ya se borraron en el generator
                    # (solo queda el directorio vacío)
                    shutil.rmtree(cleanup_dir, ignore_errors=True)
                else:
                    # Sin hit: conservar archivos para posible reintento
                    _kept_files[job_id] = {
                        'work_dir': cleanup_dir,
                        'species':  species
                    }

    threading.Thread(target=worker, daemon=True).start()
    return job_id


def get_job_queue(job_id: str):
    return _jobs.get(job_id)


def cleanup_job(job_id: str):
    _jobs.pop(job_id, None)
