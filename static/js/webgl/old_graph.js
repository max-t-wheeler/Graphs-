//create context
var context;             

function initializeContext(canvas) {

   context = canvas.getContext('experimental-webgl');

   context.viewportWidth = canvas.width;
   context.viewportHeight = canvas.height;

   if (!context) {
      alert("Could not initialize WebGL");
   }

}

//create buffers
var vertexPositions;      
var vertexPositionBuffer;

function initializeGraphVertexPositions(graph) {

   vertexPositions = [];

   for (i = 0; i < graph.order; ++i) {
      vertexPositions.push(graph.vertices[i].x);
      vertexPositions.push(graph.vertices[i].y);
      vertexPositions.push(graph.vertices[i].z);
   }

}

function initializeBuffers(graph) {

   // Create an empty buffer object to store vertex data
   vertexPositionBuffer = context.createBuffer();
   context.bindBuffer(context.ARRAY_BUFFER, vertexPositionBuffer);

   var vertexPositions = [];

   for (i = 0; i < graph.order; ++i) {
      vertexPositions.push(graph.vertices[i].x);
      vertexPositions.push(graph.vertices[i].y);
      vertexPositions.push(graph.vertices[i].z);
   }

   context.bufferData(context.ARRAY_BUFFER, new Float32Array(vertexPositions), context.STATIC_DRAW);

   vertexPositionBuffer.itemSize = graph.dimension;
   vertexPositionBuffer.numItems = graph.order;
}

//initialize matrices

var mvMatrix = mat4.create();
var pMatrix = mat4.create()

//create shaders

var vertexSrc =

   'attribute vec3 position;' + 
   'uniform mat4 Pmatrix;' +
   'uniform mat4 Vmatrix;' +
   'uniform mat4 Mmatrix;' +
   'attribute vec3 color;' +
   'varying vec3 vColor;' +

   'void main(void) {' + 
   '  gl_Position = Pmatrix*Vmatrix*Mmatrix*vec4(position, 1.0);' + 
   '  vColor = color;' +
   '}';

var fragmentSrc = 

   'precision mediump float;' +
   'varying vec3 vColor;' +

   'void main(void) {' + 
   ' gl_FragColor = vec4(vColor, 1.0);' +
   '}';

var shaderProgram;

function initializeShaders() {

   var vertexShader = context.createShader(context.VERTEX_SHADER);

   context.shaderSource(vertexShader, vertexSrc);
   context.compileShader(vertexShader);

           /////////// fragment shader ///////////

   var fragmentShader = context.createShader(context.FRAGMENT_SHADER);

   context.shaderSource(fragmentShader, fragmentSrc);
   context.compileShader(fragmentShader);

   shaderProgram = context.createProgram();

   context.attachShader(shaderProgram, vertexShader); 
   context.attachShader(shaderProgram, fragmentShader);

   context.linkProgram(shaderProgram);

   context.useProgram(shaderProgram);
}

function drawScene() {

   context.viewport(0, 0, context.viewportWidth, context.viewportHeight);
   context.clear(context.COLOR_BUFFER_BIT | context.DEPTH_BUFFER_BIT);

   mat4.perspective(45, gl.viewportWidth / gl.viewportHeight, 0.1, 100.0, pMatrix);

   mat4.identity(mvMatrix);

   context.bindBuffer(context.ARRAY_BUFFER, vertexPositionBuffer);
   context.vertexAttribPointer(shaderProgram.vertexPositionAttribute, vertexPositionBuffer.itemSize, context.FLOAT, false, 0, 0);

   setMatrixUniforms();

   context.drawArrays(context.POINTS, 0, vertexPositionBuffer.numItems);

}

function renderCanvas(element, json) {

   /*================= create canvas ==================*/

   var canvas = element;
   initializeContext(canvas);

   initializeBuffers(json);

   context.clearColor(0.0, 0.0, 0.0, 1.0);
   context.enable(context.DEPTH_TEST);

   drawScene();

}

   /*============== initialize vertices ===============*/

   /**
   var vertices = [

      -1,-1,-1, 1,-1,-1, 1, 1,-1, -1, 1,-1,
      -1,-1, 1, 1,-1, 1, 1, 1, 1, -1, 1, 1,
      -1,-1,-1, -1, 1,-1, -1, 1, 1, -1,-1, 1,
      1,-1,-1, 1, 1,-1, 1, 1, 1, 1,-1, 1,
      -1,-1,-1, -1,-1, 1, 1,-1, 1, 1,-1,-1,
      -1, 1,-1, -1, 1, 1, 1, 1, 1, 1, 1,-1, 

   ];
/**
   var colors = [

      5,3,7, 5,3,7, 5,3,7, 5,3,7,
      1,1,3, 1,1,3, 1,1,3, 1,1,3,
      0,0,1, 0,0,1, 0,0,1, 0,0,1,
      1,0,0, 1,0,0, 1,0,0, 1,0,0,
      1,1,0, 1,1,0, 1,1,0, 1,1,0,
      0,1,0, 0,1,0, 0,1,0, 0,1,0

   ];
/**
   var indices = [

      0,1,2, 0,2,3, 4,5,6, 4,6,7,
      8,9,10, 8,10,11, 12,13,14, 12,14,15,
      16,17,18, 16,18,19, 20,21,22, 20,22,23 

   ];

/**

   // Create an empty buffer object to store vertex data
   var vertexBuffer = context.createBuffer();
   context.bindBuffer(context.ARRAY_BUFFER, vertexBuffer);
   context.bufferData(context.ARRAY_BUFFER, new Float32Array(vertices), context.STATIC_DRAW);

/**
   // Create an empty buffer object to store color data
   var colorBuffer = context.createBuffer();
   context.bindBuffer(context.ARRAY_BUFFER, colorBuffer);
   context.bufferData(context.ARRAY_BUFFER, new Float32Array(colors), context.STATIC_DRAW);
/**
   // Create an empty buffer object to store index data
   var indexBuffer = context.createBuffer();
   context.bindBuffer(context.ELEMENT_ARRAY_BUFFER, indexBuffer);
   context.bufferData(context.ELEMENT_ARRAY_BUFFER, new Uint16Array(indices), context.STATIC_DRAW);

/**/

   /*============== initialize shaders ===============*/

   /**

   var vertexSrc =

      'attribute vec3 position;' + 
      'uniform mat4 Pmatrix;' +
      'uniform mat4 Vmatrix;' +
      'uniform mat4 Mmatrix;' +
      'attribute vec3 color;' +
      'varying vec3 vColor;' +

      'void main(void) {' + 
      '  gl_Position = Pmatrix*Vmatrix*Mmatrix*vec4(position, 1.0);' + 
      '  vColor = color;' +
      '}';

   var vertexShader = context.createShader(context.VERTEX_SHADER);

   context.shaderSource(vertexShader, vertexSrc);
   context.compileShader(vertexShader);

           /////////// fragment shader ///////////

   var fragmentSrc = 

   'precision mediump float;' +
   'varying vec3 vColor;' +

   'void main(void) {' + 
   ' gl_FragColor = vec4(vColor, 1.0);' +
   '}';

   var fragmentShader = context.createShader(context.FRAGMENT_SHADER);

   context.shaderSource(fragmentShader, fragmentSrc);
   context.compileShader(fragmentShader);

   var shaderProgram = context.createProgram();

   context.attachShader(shaderProgram, vertexShader); 
   context.attachShader(shaderProgram, fragmentShader);

   context.linkProgram(shaderProgram);

   /*============== associate shaders and buffers ===============*/
/**
   var Pmatrix = context.getUniformLocation(shaderProgram, "Pmatrix");
   var Vmatrix = context.getUniformLocation(shaderProgram, "Vmatrix");
   var Mmatrix = context.getUniformLocation(shaderProgram, "Mmatrix");

   // vertex buffer
   context.bindBuffer(context.ARRAY_BUFFER, vertexBuffer);

   var position = context.getAttribLocation(shaderProgram, "position");

   context.vertexAttribPointer(position, 3, context.FLOAT, false, 0, 0);
   context.enableVertexAttribArray(position);

/**

   //index buffer
   context.bindBuffer(context.ELEMENT_ARRAY_BUFFER, indexBuffer);

/**
   //color buffer
   context.bindBuffer(context.ARRAY_BUFFER, colorBuffer);

   var color = context.getAttribLocation(shaderProgram, "color");

   context.vertexAttribPointer(color, 3, context.FLOAT, false, 0, 0);
   context.enableVertexAttribArray(color);
/**
   context.useProgram(shaderProgram); 

   /*================ translation =================*/

   /**
   var Tx = 0.5, Ty = 0.5, Tz = 0.0;

   var translation = context.getUniformLocation(shader_program, 'translation');

   context.uniform4f(translation, Tx, Ty, Tz, 0.0);

   /*================== scaling ===================*/

   /**
   var Sx = 1.0, Sy = 1.5, Sz = 1.0;

   var xformMatrix = new Float32Array([
      Sx, 0.0, 0.0, 0.0,
      0.0, Sy, 0.0, 0.0,
      0.0, 0.0, Sz, 0.0, 
      0.0, 0.0, 0.0, 1.0
   ]);

   var u_xformMatrix = context.getUniformLocation(shader_program, 'u_xformMatrix');

   context.uniformMatrix4fv(u_xformMatrix, false, xformMatrix);

   /*========================= MATRIX ========================= */
/**      
   function get_projection(angle, a, zMin, zMax) {
      var ang = Math.tan((angle*0.5)*Math.PI/180);
      return [
         0.5/ang, 0 , 0, 0,
         0, 0.5*a/ang, 0, 0,
         0, 0, -(zMax+zMin)/(zMax-zMin), -1,
         0, 0, (-2*zMax*zMin)/(zMax-zMin), 0
      ];
   }

   var proj_matrix = get_projection(40, canvas.width/canvas.height, 1, 100);

   var mov_matrix = [
      1,0,0,0, 
      0,1,0,0, 
      0,0,1,0, 
      0,0,0,1
   ];

   var view_matrix = [
      1,0,0,0, 
      0,1,0,0, 
      0,0,1,0, 
      0,0,0,1
   ];

   //translating z
   view_matrix[14] = view_matrix[14]-6; //zoom

   /*================= Mouse events ======================*/
/*
   var AMORTIZATION = 0.95;
   var drag = false;
   var old_x, old_y;
   var dX = 0, dY = 0;

   var mouseDown = function(e) {
      drag = true;
      old_x = e.pageX, old_y = e.pageY;
      e.preventDefault();
      return false;
   };

   var mouseUp = function(e){
      drag = false;
   };

   var mouseMove = function(e) {
      if (!drag) return false;
      dX = (e.pageX-old_x)*2*Math.PI/canvas.width,
      dY = (e.pageY-old_y)*2*Math.PI/canvas.height;
      THETA+= dX;
      PHI+=dY;
      old_x = e.pageX, old_y = e.pageY;
      e.preventDefault();
   };

   canvas.addEventListener("mousedown", mouseDown, false);
   canvas.addEventListener("mouseup", mouseUp, false);
   canvas.addEventListener("mouseout", mouseUp, false);
   canvas.addEventListener("mousemove", mouseMove, false);

   /*=======================rotation========================*/
/*
   function rotateZ(m, angle) {
      var c = Math.cos(angle);
      var s = Math.sin(angle);
      var mv0 = m[0], mv4 = m[4], mv8 = m[8]; 
      
      m[0] = c*m[0]-s*m[1];
      m[4] = c*m[4]-s*m[5];
      m[8] = c*m[8]-s*m[9];
      m[1] = c*m[1]+s*mv0;
      m[5] = c*m[5]+s*mv4;
      m[9] = c*m[9]+s*mv8;
   }

   function rotateX(m, angle) {
      var c = Math.cos(angle);
      var s = Math.sin(angle);
      var mv1 = m[1], mv5 = m[5], mv9 = m[9];
      
      m[1] = m[1]*c-m[2]*s;
      m[5] = m[5]*c-m[6]*s;
      m[9] = m[9]*c-m[10]*s;

      m[2] = m[2]*c+mv1*s;
      m[6] = m[6]*c+mv5*s;
      m[10] = m[10]*c+mv9*s;
   }

   function rotateY(m, angle) {
      var c = Math.cos(angle);
      var s = Math.sin(angle);
      var mv0 = m[0], mv4 = m[4], mv8 = m[8];
      
      m[0] = c*m[0]+s*m[2];
      m[4] = c*m[4]+s*m[6];
      m[8] = c*m[8]+s*m[10];

      m[2] = c*m[2]-s*mv0;
      m[6] = c*m[6]-s*mv4;
      m[10] = c*m[10]-s*mv8;
   }

   /*==============graph================*/
   /*
   function render_graph(G, animation_type, vertex_type, dphi) {

      if (animation_type != "static" && animation_type != "dynamic") {

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glColor3f(1, 0, 0);

        glutSwapBuffers();

      }

       if (animation_type == "dynamic") {

           glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

           glColor3f(1, 0, 0);

           if (G.size() > 0) {

               glBegin(GL_LINES);

                   for (int i = 0; i < G.order()-1; ++i) {
                       for(int j = i+1; j < G.order(); ++j) {
                           
                           if(G.adjacency_tensor(i, j, 0) == 1) {
                               //colorNodes("white");
                               colorNodes(G.v[i]->sid+1);
                               glVertex2f(G.v[i]->x(0), G.v[i]->x(1));
                               glVertex2f(G.v[j]->x(0), G.v[j]->x(1));
                           }

                       }
                   }

               glEnd();

           }

           if (G.order() > 0) {

               glBegin(GL_POINTS);

                   for (int i = 0; i < G.order(); ++i) {

                       colorNodes(G.v[i]->sid+1);
                       
                       if (vertex_type == "dot") {
                           glVertex2f(G.v[i]->x(0), G.v[i]->x(1));
                       }
                       else if (vertex_type == "ball") {
                           filled_polygon(G.v[i]->x, G.v[i]->r, 0, 100);
                       }

                   }

               glEnd();

           }

           glFlush();

           if (dphi >= 0 ){
               dphi += sin(0.1/100)*0.1;
           }
           else {
               dphi -= 0.1/100;
           }


           glutPostRedisplay();

       }

   }
   /**/
   /*============== main ===============*/
   /**
   var THETA = 0;
   var PHI = 0;
   var time_old = 0;

   var animate = function(time) {

      var dt = time - time_old;

      if (!drag) {
         dX *= AMORTIZATION, dY*=AMORTIZATION;
         THETA+=dX, PHI+=dY;
      }
         
      mov_matrix[0] = 1, mov_matrix[1] = 0, mov_matrix[2] = 0,
      mov_matrix[3] = 0,
         
      mov_matrix[4] = 0, mov_matrix[5] = 1, mov_matrix[6] = 0,
      mov_matrix[7] = 0,
         
      mov_matrix[8] = 0, mov_matrix[9] = 0, mov_matrix[10] = 1,
      mov_matrix[11] = 0,
         
      mov_matrix[12] = 0, mov_matrix[13] = 0, mov_matrix[14] = 0,
      mov_matrix[15] = 1;

      rotateY(mov_matrix, THETA);
      rotateX(mov_matrix, PHI);
         
      time_old = time; 

      rotateZ(mov_matrix, dt*0.005);
      rotateY(mov_matrix, dt*0.002);
      rotateX(mov_matrix, dt*0.003);

      time_old = time;

      context.enable(context.DEPTH_TEST);
      context.depthFunc(context.LEQUAL);
      context.clearColor(0.5, 0.5, 0.5, 0.9);
      context.clearDepth(1.0);
      context.viewport(0.0, 0.0, canvas.width, canvas.height);
      context.clear(context.COLOR_BUFFER_BIT | context.DEPTH_BUFFER_BIT);

      context.uniformMatrix4fv(Pmatrix, false, proj_matrix);
      context.uniformMatrix4fv(Vmatrix, false, view_matrix);
      context.uniformMatrix4fv(Mmatrix, false, mov_matrix);

      context.bindBuffer(context.ELEMENT_ARRAY_BUFFER, index_buffer);
      context.drawElements(context.TRIANGLES, indices.length, context.UNSIGNED_SHORT, 0);

      window.requestAnimationFrame(animate);

   }

   animate(0);

   /**

   function draw() {

      context.enable(context.DEPTH_TEST);
      context.depthFunc(context.LEQUAL);
      context.clearColor(0.5, 0.5, 0.5, 0.9);
      context.clearDepth(1.0);
      context.viewport(0.0, 0.0, canvas.width, canvas.height);
      context.clear(context.COLOR_BUFFER_BIT | context.DEPTH_BUFFER_BIT);

   }

   /**/