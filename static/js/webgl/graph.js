var context; 

var vertexPositionBuffer;
var vertexColorBuffer;

var shaderProgram;

var vertexPositionAttribute;
var vertexColorAttribute;

var mvMatrix = mat4.create();
var pMatrix = mat4.create();
//var perspectiveMatrix;      

/******************************Quick Math***********************************/

function phi(t, n) {
   
   var p = 2*Math.PI*t/n;

   return p;

}
/*****************************Color Coding**********************************/

function colorNodes(n) {

   var mod = n%7;

   var color = [];

   switch(mod) {

      case 0:

         color = [1.0, 0.0, 0.0, 1.0];
         break;

      case 1:

         color = [0.0, 1.0, 0.0, 1.0];
         break;

      case 2:

         color = [0.0, 0.0, 1.0, 1.0];
         break;

      case 3:

         color = [1.0, 1.0, 0.0, 1.0];
         break;

      case 4:

         color = [0.0, 1.0, 1.0, 1.0];
         break;

      case 5:

         color = [1.0, 0.0, 1.0, 1.0];
         break;

      case 6:

         color = [1.0, 1.0, 1.0, 1.0];
         break;

   }

    return color;

}   

/********************************Shapes*************************************/

function generatePolygonVertices(vertices, r, radius, numNodes) {

   for (i = 1; i <= numNodes; ++i) {
      
      var x = radius*Math.sin(phi(i, numNodes)) + r[0];
      var y = radius*Math.cos(phi(i, num_nodes))+ r[1];
      var z = 0.0;
      
      vertices.push(x);
      vertices.push(y);
      vertices.push(z);
   
   }

}

/***************************Matrix Utilities********************************/
/*
function loadIdentity() {
   mvMatrix = Matrix.I(4);
}

function multMatrix(m) {
   mvMatrix = mvMatrix.x(m);
}

function mvTranslate(v) {
   multMatrix(Matrix.Translation($V([v[0], v[1], v[2]])).ensure4x4());
}

function setMatrixUniforms() {

   var pUniform = context.getUniformLocation(shaderProgram, "uPMatrix");
   context.uniformMatrix4fv(pUniform, false, new Float32Array(pMatrix.flatten()));

   var mvUniform = context.getUniformLocation(shaderProgram, "uMVMatrix");
   context.uniformMatrix4fv(mvUniform, false, new Float32Array(mvMatrix.flatten()));

}
*/

function setMatrixUniforms() {

   context.uniformMatrix4fv(shaderProgram.pMatrixUniform, false, pMatrix);
   context.uniformMatrix4fv(shaderProgram.mvMatrixUniform, false, mvMatrix);

}

/************************WebGL***********************************/

function initializeContext(canvas) {

   context  = null;

   try {
      
      context = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
      
      context.viewportWidth = canvas.width;
      context.viewportHeight = canvas.height;
   
   }
   catch(e) {

   }

   if (!context) {
      alert("Could not initialize WebGL");
   }

}

function getShader(context, id, type) {

   var shaderScript, shaderSource, currentChild, shader;

   shaderScript = document.getElementById(id);

   if (!shaderScript) {
      console.log('test');
      return null;
   }

   shaderSource = shaderScript.text;

   if (!type) {

      if (shaderScript.type == 'x-shader/x-fragment') {
         type = context.FRAGMENT_SHADER;
      }
      else if ('x-shader/x-vertex') {
         type = context.VERTEX_SHADER;
      }
      else {
         return null;
      }

   }

   shader = context.createShader(type);

   context.shaderSource(shader, shaderSource);
   context.compileShader(shader);

   if(!context.getShaderParameter(shader, context.COMPILE_STATUS)) {

      console.log('An error occurred compiling the shaders: ' + context.getShaderInfoLog(shader));
      
      context.deleteShader(shader);

      return null;
   
   }

   return shader;

}

function initializeShaders() {

   var fragmentShader = getShader(context, 'shader-fs');
   var vertexShader = getShader(context, 'shader-vs');

   shaderProgram = context.createProgram();

   context.attachShader(shaderProgram, vertexShader); 
   context.attachShader(shaderProgram, fragmentShader);
   context.linkProgram(shaderProgram);

   if (!context.getProgramParameter(shaderProgram, context.LINK_STATUS)) {
      console.log('Unable to initialize the shader program: ' + context.getProgramInfoLog(shaderProgram));
   }

   context.useProgram(shaderProgram);

   vertexPositionAttribute = context.getAttribLocation(shaderProgram, 'aVertexPosition');
   context.enableVertexAttribArray(vertexPositionAttribute);

   vertexColorAttribute = context.getAttribLocation(shaderProgram, 'aVertexColor');
   context.enableVertexAttribArray(vertexColorAttribute);

   shaderProgram.pMatrixUniform = context.getUniformLocation(shaderProgram, "uPMatrix");
   shaderProgram.mvMatrixUniform = context.getUniformLocation(shaderProgram, "uMVMatrix");

}

function initializeBuffers(graph) {

   var dimension = graph.dimension;
   var numNodes = 100;
   var numVertices = graph.order;

   // Create an empty buffer object to store vertex data
   vertexPositionBuffer = context.createBuffer();

   var vertexPositions = [];
   var vertexColors = [];

/**
   vertexPositions = [
      1.0,  1.0,  0.0,
      -1.0, 1.0,  0.0,
      1.0,  -1.0, 0.0,
      -1.0, -1.0, 0.0      
   ];

/**

   for (i = 0; i < graph.order; ++i) {

      vertexPositions.push(graph.vertices[i].x);
      vertexPositions.push(graph.vertices[i].y);
      vertexPositions.push(graph.vertices[i].z);
   
   }

/**/

   for (i = 0; i < numVertices; ++i) {
      
      var vertexRadius = graph.vertices[i].r;

      for (j = 1; j <= numNodes; ++j) {

         //var x = graph.vertices[i].x;
         //var y = graph.vertices[i].y;

         var x = vertexRadius*Math.sin(phi(j, numNodes)) + graph.vertices[i].x;
         var y = vertexRadius*Math.cos(phi(j, numNodes)) + graph.vertices[i].y;
         var z = graph.vertices[i].z;

         vertexPositions.push(x);
         vertexPositions.push(y);
         vertexPositions.push(z);
      
         vertexColors = vertexColors.concat(colorNodes(graph.vertices[i].sid));

      }

   }

/**/

   context.bindBuffer(context.ARRAY_BUFFER, vertexPositionBuffer);
   context.bufferData(context.ARRAY_BUFFER, new Float32Array(vertexPositions), context.STATIC_DRAW);

   vertexPositionBuffer.itemSize = dimension;
   vertexPositionBuffer.numItems = numVertices*numNodes;

/**

   vertexColors = [
    1.0,  1.0,  1.0,  1.0,    // white
    1.0,  0.0,  0.0,  1.0,    // red
    0.0,  1.0,  0.0,  1.0,    // green
    0.0,  0.0,  1.0,  1.0     // blue
  ];

/**

   for (i = 0; i < graph.order; ++i) {
      vertexColors = vertexColors.concat(colorNodes(graph.vertices[i].sid));
   }

/**/

  vertexColorBuffer = context.createBuffer();
  context.bindBuffer(context.ARRAY_BUFFER, vertexColorBuffer);
  context.bufferData(context.ARRAY_BUFFER, new Float32Array(vertexColors), context.STATIC_DRAW);

}

function drawScene() {

   context.viewport(0, 0, context.viewportWidth, context.viewportHeight);
   context.clear(context.COLOR_BUFFER_BIT | context.DEPTH_BUFFER_BIT);

/*
   perspectiveMatrix = makePerspective(45, 1.0, 0.1, 100.0);

   loadIdentity();
   mvTranslate([-0.0, 0.0, -6.0]);
*/

   mat4.perspective(pMatrix, 45, context.viewportWidth/context.viewportHeight, 0.1, 100);
   mat4.identity(mvMatrix);

   mat4.translate(mvMatrix, mvMatrix, [-0.0, 0.0, -6.0]);

   context.bindBuffer(context.ARRAY_BUFFER, vertexPositionBuffer);
   context.vertexAttribPointer(vertexPositionAttribute, 3, context.FLOAT, false, 0, 0);

   context.bindBuffer(context.ARRAY_BUFFER, vertexColorBuffer);
   context.vertexAttribPointer(vertexColorAttribute, 4, context.FLOAT, false, 0, 0);
   
   setMatrixUniforms();

   context.drawArrays(context.POINTS, 0, vertexPositionBuffer.numItems);
   //context.drawArrays(context.TRIANGLE_STRIP, 0, 4);

}

function tick() {

   requestAnimFrame(tick);
   drawScene();
   //animate();

}

function renderCanvas(element, json) {

   var canvas = element;

   initializeContext(canvas);

   if (!context) {
      return;
   }

   initializeShaders();
   initializeBuffers(json);

   context.clearColor(0.0, 0.0, 0.0, 1.0);
   context.enable(context.DEPTH_TEST);
   context.depthFunc(context.LEQUAL);

   tick();

} 