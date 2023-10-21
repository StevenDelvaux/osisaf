from waitress import serve
import osisaf-api
serve(osisaf-api.app, host='0.0.0.0', port=8080)