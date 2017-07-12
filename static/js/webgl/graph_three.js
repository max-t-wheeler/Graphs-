var camera, scene, renderer;

var geometry, group; 

var dimension, numVertices, numEdges, radius;

function init(graph) {
   
   var container = document.createElement('div');

   container.id = 'graph-div';

   document.body.appendChild(container);

   camera = new THREE.PerspectiveCamera(60, 1, 1, 1000);

   camera.position.z = 10;

   scene = new THREE.Scene();

   dimension = graph.dimension;
   numVertices = graph.order;
   numEdges = graph.size;
   radius = graph.radius;

   if (numEdges > 0) {

      for (var i = 0; i < numVertices - 1; ++i) {

         for (var j = i + 1; j < numVertices; ++j) {

            if (graph.adjacency_tensor[i].degree_sequence[j] == 1) {
               
               var material = new THREE.LineBasicMaterial({ color: 0xffffff });
               var geometry = new THREE.Geometry();

               geometry.vertices.push(
                  new THREE.Vector3(graph.vertices[i].x, graph.vertices[i].y, graph.vertices[i].z), 
                  new THREE.Vector3(graph.vertices[j].x, graph.vertices[j].y, graph.vertices[j].z) 
               );

               var line = new THREE.Line(geometry, material);

               scene.add(line);

            }

         }

      }

   }

   if (numVertices > 0) {

      geometry = new THREE.Geometry();
      group = new THREE.Group();

      for (var i = 0; i < numVertices; ++i) {

   /**/
         var material = new THREE.SpriteCanvasMaterial( { 
            
            color: colorNodes(graph.vertices[i].sid),
            
            program: function(context) {

               context.beginPath();
               context.arc(0, 0, graph.vertices[i].r, 0, 2 * Math.PI, true);
               context.fill();

            }  
            

         });

         var vertex = new THREE.Sprite(material);

         vertex.position.x = graph.vertices[i].x;
         vertex.position.y = graph.vertices[i].y;

         if (dimension == 3) {
            vertex.position.z = graph.vertices[i].z;
         }

         vertex.position.normalize();
         vertex.scale.multiplyScalar(2);

         group.add(vertex);


   /**
         geometry = new THREE.SphereGeometry(graph.vertices[i].r, 32, 32);
         var material = new THREE.MeshBasicMaterial({ color: colorNodes(graph.vertices[i].sid) });

         var mesh = new THREE.Mesh(geometry, material);

         mesh.position.x = graph.vertices[i].x;
         mesh.position.y = graph.vertices[i].y;

         if (dimension == 3) {
            mesh.position.z = graph.vertices[i].z;
         }

         group.add(mesh);
   /**/

      }

      scene.add(group);

   }


   renderer = new THREE.WebGLRenderer();
   renderer.setClearColor( 0x000000 );
   renderer.setPixelRatio( window.devicePixelRatio );
   renderer.setSize(450, 450);
   container.appendChild( renderer.domElement );

   var controls = new THREE.OrbitControls(camera, renderer.domElement);
   controls.target.set(0, 0, 0);
   controls.update();


}

function render() {
   renderer.render(scene, camera);
}

var dt = 0;


function animate() {
   
   for (var i = 0; i < group.children.length; ++i) {
      group.children[i].position.x += group.children[i].position.x*Math.sin(i*dt/2)/50;
      group.children[i].position.y -= group.children[i].position.y*Math.sin(i*dt/2)/50;
   }

   dt += 0.01;

   if (dt > 2*Math.PI*1e6) {
      dt -= 2*Math.PI*1e6;
   }

}

function tick() {

   requestAnimationFrame(tick);
   //animate();
   render();
}