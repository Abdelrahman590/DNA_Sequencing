// Chart Configuration
const chartConfig = {
    options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
            legend: {
                position: 'top',
            },
            title: {
                display: true,
                text: 'تحليل التسلسل'
            }
        }
    }
};

// Create Nucleotide Composition Chart
function createCompositionChart(data) {
    const ctx = document.getElementById('compositionChart');
    if (!ctx) return;

    new Chart(ctx, {
        type: 'doughnut',
        data: {
            labels: ['A', 'T', 'G', 'C'],
            datasets: [{
                data: [
                    data.composition.A,
                    data.composition.T,
                    data.composition.G,
                    data.composition.C
                ],
                backgroundColor: [
                    '#FF6384',
                    '#36A2EB',
                    '#FFCE56',
                    '#4BC0C0'
                ]
            }]
        },
        options: {
            ...chartConfig.options,
            plugins: {
                ...chartConfig.options.plugins,
                title: {
                    ...chartConfig.options.plugins.title,
                    text: 'تركيب النيوكليوتيدات'
                }
            }
        }
    });
}

// Create GC Content Chart
function createGCChart(gcContent) {
    const ctx = document.getElementById('gcChart');
    if (!ctx) return;

    new Chart(ctx, {
        type: 'gauge',
        data: {
            datasets: [{
                value: gcContent,
                minValue: 0,
                maxValue: 100,
                backgroundColor: ['#FF6384', '#36A2EB']
            }]
        },
        options: {
            ...chartConfig.options,
            plugins: {
                ...chartConfig.options.plugins,
                title: {
                    ...chartConfig.options.plugins.title,
                    text: 'محتوى GC'
                }
            }
        }
    });
}

// Create Statistics Chart
function createStatsChart(data) {
    const ctx = document.getElementById('statsChart');
    if (!ctx) return;

    new Chart(ctx, {
        type: 'bar',
        data: {
            labels: ['AT Content', 'GC Content', 'Purine Content', 'Pyrimidine Content'],
            datasets: [{
                label: 'النسبة المئوية',
                data: [
                    data.statistics.at_content,
                    data.gc_content,
                    data.statistics.purine_content,
                    data.statistics.pyrimidine_content
                ],
                backgroundColor: [
                    '#FF6384',
                    '#36A2EB',
                    '#FFCE56',
                    '#4BC0C0'
                ]
            }]
        },
        options: {
            ...chartConfig.options,
            scales: {
                y: {
                    beginAtZero: true,
                    max: 100
                }
            }
        }
    });
}

// File Upload Handling
document.addEventListener('DOMContentLoaded', function() {
    const fileInput = document.getElementById('file-input');
    const uploadArea = document.querySelector('.upload-area');
    
    if (fileInput && uploadArea) {
        // Drag and drop functionality
        uploadArea.addEventListener('dragover', (e) => {
            e.preventDefault();
            uploadArea.classList.add('dragover');
        });

        uploadArea.addEventListener('dragleave', () => {
            uploadArea.classList.remove('dragover');
        });

        uploadArea.addEventListener('drop', (e) => {
            e.preventDefault();
            uploadArea.classList.remove('dragover');
            
            if (e.dataTransfer.files.length) {
                fileInput.files = e.dataTransfer.files;
                handleFileSelection();
            }
        });

        // Regular file input handling
        fileInput.addEventListener('change', handleFileSelection);
    }
});

// Handle file selection
function handleFileSelection() {
    const fileInput = document.getElementById('file-input');
    const fileName = document.getElementById('file-name');
    const uploadForm = document.getElementById('upload-form');
    
    if (fileInput.files.length) {
        fileName.textContent = fileInput.files[0].name;
        uploadForm.classList.add('has-file');
    }
}

// Initialize charts when results are available
function initializeCharts(data) {
    createCompositionChart(data);
    createGCChart(data.gc_content);
    createStatsChart(data);
}

// Export results to various formats
function exportResults(format) {
    const resultsDiv = document.getElementById('analysis-results');
    if (!resultsDiv) return;

    switch(format) {
        case 'pdf':
            // Add PDF export logic
            break;
        case 'csv':
            // Add CSV export logic
            break;
    }
}
