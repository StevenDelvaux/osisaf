from waitress import serve
import osisaf
serve(osisaf.app, host='0.0.0.0', port=8080)