$(document).ready(function(){

    $('#set-num-vert').click(function() {

        var entry = $('#num-species').val();

        $('#num-vert').empty();

        if (entry >= 10) {
           alert("Number of species must be less than to 10");
        }
        else {

            for (i = 1; i <= entry; ++i) {

               var input = document.createElement('input');

               var label = document.createElement('Label');
                
               label.innerHTML = '<br> Number of vertices of Species ' + i + '<br>';

               input.setAttribute('type', 'number');
               input.setAttribute('name', 'num-vert-' + i);
               input.setAttribute('class', 'num-vert')
               input.setAttribute('style', 'width:200px');

               label.setAttribute('style', 'font-weight:normal');

               var num_vert = document.getElementById('num-vert');

               num_vert.appendChild(label);
               num_vert.appendChild(input);

            }

        }

    });

    $('#graph-dimension-2D').click(function() {

        var selection = document.getElementsByName('coords-2')
        
        $(selection).remove();

    });

    $('#graph-dimension-3D').click(function() {
        
        var radio = document.getElementsByName('graph-center');

        if (radio[1].checked) {
            
            $('#coords').empty();

            for (i = 0; i < 3; ++i) {

                var input = document.createElement('input');

                var label = document.createElement('Label');
                label.innerHTML = '<br>';


                input.setAttribute('type', 'number');
                input.setAttribute('name', 'coords-' + i);
                input.setAttribute('style', 'width:200px');

                var coords = document.getElementById('coords');

                coords.appendChild(label);
                coords.appendChild(input);

            }      

        }

    });

    $('#center-origin').click(function() {

        $('#coords').empty();

    });

    $('#center-custom').click(function() {

        $('#coords').empty();

        var radio = document.getElementsByName('graph-dimension');

        var dimension;

        for (i = 0; i < radio.length; ++i) {
            if (radio[i].checked) {
                dimension = radio[i].value;
            }
        }

        for (i = 0; i < dimension; ++i) {

            var input = document.createElement('input');

            var label = document.createElement('Label');
            label.innerHTML = '<br>';


            input.setAttribute('type', 'number');
            input.setAttribute('name', 'coords-' + i);
            input.setAttribute('style', 'width:200px');

            var coords = document.getElementById('coords');

            coords.appendChild(label);
            coords.appendChild(input);

        }

    });

/**
    $('#graph-inputs').submit(function(event) {
        event.preventDefault();
        event.stopPropagation();
        fetch('/graph_inputs', {
            method: 'post',
            body: new FormData(document.getElementById('graph-inputs'))
        }).then(function(response) {
            return response.json();
        }).then(function(graph){
            console.log(graph);

            initialize_webgl();

            /**

            initialize_context();           

            var vertex_positions = [];

            for (i = 0; i < graph.order; ++i) {
                vertex_positions.push(graph.vertices[i].x);
                vertex_positions.push(graph.vertices[i].y);
                vertex_positions.push(graph.vertices[i].z);
            }

            context.enable(context.DEPTH_TEST);
            context.depthFunc(context.LEQUAL);
            context.clearColor(0.5, 0.5, 0.5, 0.9);
            context.clearDepth(1.0);
            context.viewport(0.0, 0.0, canvas.width, canvas.height);
            context.clear(context.COLOR_BUFFER_BIT | context.DEPTH_BUFFER_BIT);

            /**

        });
        return false;
    });

/**/

    $('#reset').click(function() {
        $('#num-vert').empty();
        $('#coords').empty();
    });

});