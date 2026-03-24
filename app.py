import os
import io
import re
import base64
from collections import Counter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from flask import Flask, render_template, request, jsonify

app = Flask(__name__)


def parse_numbers(text):
    """Parse numbers from text with any common separator (comma, semicolon, space, tab, newline)."""
    parts = re.split(r'[,;\s\t\n\r]+', text.strip())
    numbers = []
    for part in parts:
        part = part.strip()
        if part:
            try:
                numbers.append(float(part))
            except ValueError:
                raise ValueError(f"'{part}' no es un número válido")
    return numbers


def process_survival_data(data):
    """Process survival data and generate KM curve."""
    counter = Counter(data)
    sorted_data = [(count, value) for value, count in counter.items()]
    sorted_data.sort(key=lambda x: x[1])

    total_animals = sum(count for count, _ in sorted_data)
    animals_at_risk = total_animals
    survival_probability = 1.0
    times = [0]
    survivals = [1.0]
    table_rows = []

    for deaths, time in sorted_data:
        survival_probability *= (animals_at_risk - deaths) / animals_at_risk
        times.extend([time, time])
        survivals.extend([survivals[-1], survival_probability])
        table_rows.append({
            'deaths': int(deaths),
            'time': time,
            'at_risk': int(animals_at_risk),
            'survival': round(survival_probability, 4)
        })
        animals_at_risk -= deaths

    plt.figure(figsize=(10, 6))
    plt.plot(times, survivals, 'b-', linewidth=2)
    plt.scatter(times[1::2], survivals[1::2], color='blue', zorder=5)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlabel('Tiempo (días)', fontsize=12)
    plt.ylabel('Probabilidad de supervivencia', fontsize=12)
    plt.title('Curva de Kaplan-Meier', fontsize=14)
    plt.ylim(-0.05, 1.05)

    img_bytes = io.BytesIO()
    plt.savefig(img_bytes, format='png', dpi=300, bbox_inches='tight')
    plt.close()
    img_bytes.seek(0)

    text_output = "Muertes\tTiempo\tAnimales en riesgo\tProbabilidad de supervivencia\n"
    for row in table_rows:
        text_output += f"{row['deaths']}\t{row['time']}\t{row['at_risk']}\t{row['survival']:.4f}\n"

    return img_bytes, text_output, table_rows


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/kaplan-meier', methods=['GET', 'POST'])
def kaplan_meier():
    if request.method == 'GET':
        return render_template('kaplan_meier.html')

    try:
        numbers_text = ''
        if 'file' in request.files and request.files['file'].filename:
            file = request.files['file']
            numbers_text = file.read().decode('utf-8')
        else:
            numbers_text = request.form.get('numbers', '')

        if not numbers_text.strip():
            return jsonify({'error': 'No se han introducido datos. Por favor, introduce los tiempos de supervivencia.'}), 400

        numbers = parse_numbers(numbers_text)
        if not numbers:
            return jsonify({'error': 'No se encontraron números válidos en los datos introducidos.'}), 400

        img_bytes, text_output, table_rows = process_survival_data(numbers)
        img_b64 = base64.b64encode(img_bytes.getvalue()).decode('utf-8')

        return jsonify({
            'image': img_b64,
            'text_output': text_output,
            'table': table_rows
        })

    except ValueError as e:
        return jsonify({'error': str(e)}), 400
    except Exception as e:
        return jsonify({'error': f'Error inesperado: {str(e)}'}), 500


if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=False)
