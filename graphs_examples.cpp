#include <time.h>

#include "OpenGL_Utils.h"

unsigned int n = 1;
float dphi = 0.01;

arma::uvec g_vertex_set_sizes = {200, 200, 200, 200, 200, 200};
arma::uvec g_vec_vertex_set_sizes = {5, 5, 5, 5, 5, 5};
arma::vec center = arma::zeros(3);

graph g(0, g_vertex_set_sizes, 1, center, 100, 1);

std::vector<graph*> g_vec {};


void Graph() {

/**/
    g.initialize("null", "partite");

    for (int i = 0; i < g.v.size(); ++i) {

        //vertex set 3 [null, concentric]
        //g.v[i]->x(0) -= g.v[i]->x(0)*sin(i*dphi/2);
        //g.v[i]->x(1) += g.v[i]->x(1)*sin(i*dphi/2);

        //vertex set 3 (absurd) [r=0.0001, null, concentric] [additive]
        g.v[i]->x(0) += cos(10*dphi)*g.v[i]->x(0)*sin(i*dphi/2);
        g.v[i]->x(1) += cos(10*dphi)*g.v[i]->x(1)*sin(i*dphi/2);

        //vertex set 3 (absurd) [r=0.0001, null, concentric] [subtractive]
        //g.v[i]->x(0) -= cos(10*dphi)*g.v[i]->x(0)*sin(i*dphi/2);
        //g.v[i]->x(1) += cos(10*dphi)*g.v[i]->x(1)*sin(i*dphi/2);
    }

   render_graph(g, "dynamic", "ball", dphi);

}

void Graphs() {

    for (int i = 0; i < n; ++i){

        double r = 1.0; 

        if (i%2==0) {
            r = 1-(5*sin(10*dphi));
        }
        else {
            r = 1 - (5*cos(10*dphi));
        }

        graph * g = new graph(i, g_vec_vertex_set_sizes, r, center, 100, 1);;     //create new financial intermediary
        
        g->initialize("complete k-partite", "concentric");  
        
        g_vec.push_back(g);                      
    }

    render_graphs(g_vec, "dynamic", "dot", dphi);

    for (int i = 0; i < g_vec.size(); ++i){
            delete g_vec[i];
        }

    g_vec.clear(); 

}

void init(void) {

    glClearColor(0.0, 0.0, 0.0, 0.0);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

}

int main(int argc, char ** argv)
{       

    std::string query;

    output("Enter graph or graphs");

    std::cin >> query;

    srand(time(NULL));
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowPosition(400, 10);
    glutInitWindowSize(900,900);
    glutCreateWindow("Example");

    init();

    if (query == "graph") {
        glutDisplayFunc(Graph);
    }
    else if (query == "graphs") {
        glutDisplayFunc(Graphs);
    }
    else {

        output("Entry specified is invalid");
        
        return 0;
    
    }

    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);

    glutMainLoop();
    
    return 0;
}
