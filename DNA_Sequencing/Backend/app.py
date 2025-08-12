from flask import Flask, render_template, request, jsonify, flash, redirect, url_for, Response
import io
import csv
from reportlab.pdfgen import canvas
import os
from werkzeug.utils import secure_filename
from analyzer import DNAAnalyzer
import json
from datetime import datetime
import logging
import re
import secrets

# إعداد logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__, 
                template_folder='../Frontend/Templets',
                static_folder='../static')
app.secret_key = secrets.token_hex(16)  # مفتاح آمن
app.config['UPLOAD_FOLDER'] = os.path.join(os.path.dirname(__file__), 'static/uploads')
app.config['RESULTS_FOLDER'] = 'static/results'
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16MB max file size

# إنشاء المجلدات المطلوبة
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['RESULTS_FOLDER'], exist_ok=True)

# الصيغ المسموحة للملفات
ALLOWED_EXTENSIONS = {'fasta', 'fa', 'fna', 'fastq', 'fq', 'txt', 'seq', 'ffn', 'faa'}

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

def validate_dna_sequence(sequence):
    """التحقق من صحة تسلسل DNA"""
    if len(sequence) < 10:
        return False, "التسلسل قصير جداً (الحد الأدنى 10 نيوكليوتيد)"
    
    if not re.match(r'^[ATGCN]+$', sequence.upper()):
        return False, "التسلسل يحتوي على رموز غير صالحة (يُسمح فقط بـ A, T, G, C, N)"
    
    return True, ""

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        analyzer = DNAAnalyzer()
        
        try:
            # فحص إذا كان هناك ملف مرفوع
            if 'file' in request.files:
                file = request.files['file']
                if file.filename != '' and allowed_file(file.filename):
                    # إنشاء اسم ملف فريد
                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    filename = f"{timestamp}_{secure_filename(file.filename)}"
                    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                    
                    file.save(filepath)
                    logger.info(f"File uploaded: {filename}")
                    
                    results = analyzer.analyze_from_file(filepath)
                    return render_template('result.html', results=results, filename=filename)
                else:
                    flash('يرجى اختيار ملف بصيغة صحيحة (FASTA, FASTQ, TXT)')
                    return redirect(url_for('index'))
            
            # أو إذا كان النص مُدخل مباشرة
            elif 'sequence_text' in request.form:
                sequence = request.form['sequence_text'].strip().upper()
                
                # التحقق من صحة التسلسل
                is_valid, error_msg = validate_dna_sequence(sequence)
                if not is_valid:
                    flash(error_msg)
                    return redirect(url_for('index'))
                
                results = analyzer.analyze_sequence(sequence)
                return render_template('result.html', results=results, filename="إدخال مباشر")
            
            else:
                flash('يرجى اختيار ملف أو إدخال تسلسل DNA')
                return redirect(url_for('index'))
                
        except Exception as e:
            logger.error(f"Error in upload_file: {str(e)}")
            flash(f'حدث خطأ أثناء المعالجة: {str(e)}')
            return redirect(url_for('index'))
    
    # GET request - redirect to index
    return redirect(url_for('index'))

@app.route('/api/analyze', methods=['POST'])
def api_analyze():
    """API endpoint للتحليل السريع"""
    try:
        if not request.is_json:
            return jsonify({'error': 'Content-Type يجب أن يكون application/json'}), 400
        
        data = request.get_json()
        sequence = data.get('sequence', '').strip().upper()
        
        # التحقق من صحة التسلسل
        is_valid, error_msg = validate_dna_sequence(sequence)
        if not is_valid:
            return jsonify({'error': error_msg}), 400
        
        analyzer = DNAAnalyzer()
        results = analyzer.analyze_sequence(sequence)
        
        logger.info(f"API analysis completed for sequence length: {len(sequence)}")
        return jsonify(results)
    
    except Exception as e:
        logger.error(f"Error in API analysis: {str(e)}")
        return jsonify({'error': f'حدث خطأ أثناء التحليل: {str(e)}'}), 500

# Headers أمان
@app.after_request
def after_request(response):
    response.headers['X-Content-Type-Options'] = 'nosniff'
    response.headers['X-Frame-Options'] = 'DENY'
    response.headers['X-XSS-Protection'] = '1; mode=block'
    return response

@app.route('/export/<format>', methods=['POST'])
def export_results(format):
    """تصدير النتائج بصيغ مختلفة"""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No data provided'}), 400

        if format == 'csv':
            # Generate CSV
            output = io.StringIO()
            writer = csv.writer(output)
            
            # Write headers
            writer.writerow(['Property', 'Value'])
            
            # Basic sequence info
            writer.writerow(['Sequence Length', data['sequence_info']['length']])
            writer.writerow(['GC Content', f"{data['gc_content']}%"])
            
            # Composition
            for base, percentage in data['composition'].items():
                writer.writerow([f'{base} Content', f"{percentage}%"])
            
            # Statistics
            for stat, value in data['statistics'].items():
                if stat != 'summary':
                    writer.writerow([stat.replace('_', ' ').title(), value])
            
            return Response(
                output.getvalue(),
                mimetype='text/csv',
                headers={'Content-Disposition': 'attachment; filename=sequence_analysis.csv'}
            )
            
        elif format == 'pdf':
            # Create PDF
            buffer = io.BytesIO()
            pdf = canvas.Canvas(buffer)
            
            # Title
            pdf.setFont("Helvetica-Bold", 16)
            pdf.drawString(50, 800, "DNA Sequence Analysis Report")
            
            # Basic Info
            pdf.setFont("Helvetica", 12)
            y = 750
            pdf.drawString(50, y, f"Sequence Length: {data['sequence_info']['length']}")
            y -= 20
            pdf.drawString(50, y, f"GC Content: {data['gc_content']}%")
            
            # Composition
            y -= 40
            pdf.setFont("Helvetica-Bold", 14)
            pdf.drawString(50, y, "Nucleotide Composition")
            pdf.setFont("Helvetica", 12)
            for base, percentage in data['composition'].items():
                y -= 20
                pdf.drawString(70, y, f"{base}: {percentage}%")
            
            # Statistics
            y -= 40
            pdf.setFont("Helvetica-Bold", 14)
            pdf.drawString(50, y, "Statistics")
            pdf.setFont("Helvetica", 12)
            for stat, value in data['statistics'].items():
                if stat != 'summary':
                    y -= 20
                    pdf.drawString(70, y, f"{stat.replace('_', ' ').title()}: {value}")
            
            pdf.save()
            buffer.seek(0)
            
            return Response(
                buffer.getvalue(),
                mimetype='application/pdf',
                headers={'Content-Disposition': 'attachment; filename=sequence_analysis.pdf'}
            )
            
        else:
            return jsonify({'error': 'Unsupported format'}), 400
            
    except Exception as e:
        logger.error(f"Error in export: {str(e)}")
        return jsonify({'error': str(e)}), 500

# Error handlers
@app.errorhandler(413)
def too_large(error):
    flash('الملف كبير جداً. الحد الأقصى 16MB')
    return redirect(url_for('upload_file'))

@app.errorhandler(404)
def not_found_error(error):
    return render_template('index.html'), 404  # إعادة توجيه للرئيسية

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
