from flask import Flask, render_template, request, jsonify
import subprocess
import os

app = Flask(__name__)
ret = {}

@app.route('/')
def index():
   return render_template('index.html')

@app.route('/graph_inputs', methods=['POST'])
def processInput():

    num_species = request.form.get('num-species')

    ret = []
    #ret = {}
	
    ret.append(num_species)
    #ret['n_sp'] = num_species

    #print('num_species')
    #print(num_species)

    int_num_species = int(num_species)

    num_vert = ''

    #print('num_vert')

    for i in range(0, int_num_species):
    	num_vert = request.form.get('num-vert-' + chr(i+49)) 
        ret.append(num_vert) 
    	#print(num_vert)

    graph_dimension = request.form.get('graph-dimension')

    ret.append(graph_dimension)

    #print('graph_dimension')
    #print(graph_dimension)

    int_graph_dimension = int(graph_dimension)

    graph_center = request.form.get('graph-center')

    ret.append(graph_center)

    #print('graph_center')
    #print(graph_center)

    coords = ''
    #print('coords')

    for i in range(0, int_graph_dimension):
    	
        if graph_center == 'custom':
            coords = coords + request.form.get('coords-' + chr(i+48))
            ret.append(coords[i]) 
            #print(coords[i])

        #else:
        #coords.append(chr(48))
        #ret.append(coords[i]) 
        #print(coords[i])

    #print(coords)


    graph_radius = request.form.get('graph-radius')
    graph_layout = request.form.get('graph-layout')
    graph_type = request.form.get('graph-type')

    ret.append(graph_radius)
    ret.append(graph_layout)
    ret.append(graph_type)

    #print('graph_radius')
    #print(graph_radius)
    #print('graph_layout')
    #print(graph_layout)
    #print('graph_type')
    #print(graph_type)

    #print(ret)

    #print(jsonify(ret))

    output = ''

    for i in range(0, len(ret)-1):
        output = output + ret[i] + '+'
    
    output = output + ret[len(ret)-1]

    
    #print(output)

    call = subprocess.check_output(
    	[
	    	#test
	    	"./generator",
            output 
	    	#num_species,
	    	#num_vert,
	    	#graph_dimension,
	    	#graph_center,
	    	#graph_radius,
	    	#graph_layout,
	    	#graph_type
    	]
    )

    #print(call)

    return call


@app.route('/simulation', methods = ['POST', 'GET'])
def simulation():

	if request.method == 'POST':
		result = request.form 
		return render_template('simulation.html')

if __name__ == '__main__':
	app.run(debug = True)