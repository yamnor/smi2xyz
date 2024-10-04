from bottle import Bottle, template

app = Bottle()

@app.route('/')
def index():
  return '<b>Hello World</b>!'

@app.route('/hello/<name>')
def index(name):
  return template('<b>Hello {{name}}</b>!', name=name)
