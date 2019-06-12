#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <windows.h>
#include <iostream>
#include <fstream>
using namespace std;
#include<mpi.h>

#include<GLFW/glfw3.h>

//#pragma comment(lib,"glew32.lib")  
//#pragma comment(lib,"glut32.lib")  
#pragma comment(lib,"OpenGL32.lib")
#pragma comment(lib,"glfw3.dll")

//Macros to use 1D arrays with ease
#define PX(i) (3*i+1)
#define PY(i) (3*i+2)
#define MASS(i) (3*i+3)

#define VX(i) (4*i+0)
#define VY(i) (4*i+1)
#define AX(i) (4*i+2)
#define AY(i) (4*i+3)

//Gravitational constant
double G=0.0001;
//Fixed delta t
double dt=0.005;
//Barnes Hut critical distance
double rcutoff=0.35;

//To avoid almost infinity values in force calculation
double rlimit=0.03;

//Nodes of the Barnes Hut Quadtree
struct Node{
    struct Node *children[4];
    int external;
	
	//Center of Mass cordinates
    double CMX;
    double CMY;
    double mass;

	//Top right corner coordinates
    double TRX;
    double TRY;

	//Lower left corner coordinates
    double LLX;
    double LLY;

	//Geometrical center, tradeoff of space for speed
    double GCX;
    double GCY;
};

void buildTree(struct Node* node, double* shrdBuff, int *indexes, int n){
    if(n==1){ //This is an external node!
        node->external=1;
        node->CMX=shrdBuff[PX(indexes[0])];
        node->CMY=shrdBuff[PY(indexes[0])];
        node->mass=shrdBuff[MASS(indexes[0])];
    } else {
        node->external=0;
		//Arrays that will hold the indexes of the particles in each node
        int *NEi = (int *) malloc(sizeof(int)*n);
        int *NWi = (int *) malloc(sizeof(int)*n);
        int *SWi = (int *) malloc(sizeof(int)*n);
        int *SEi = (int *) malloc(sizeof(int)*n);
		//Counter for those particles
        int NWc=0, SWc=0,SEc=0, NEc=0;

        int i;
        for(i=0;i<n;i++){
            if(shrdBuff[PY(indexes[i])] < node->GCY ){ //South moitie
                if(shrdBuff[PX(indexes[i])] < node->GCX){ //West wing
                    SWi[SWc]=indexes[i];
                    SWc++;
                } else {
                    SEi[SEc]=indexes[i];
                    SEc++;
                }
            } else {
                if(shrdBuff[PX(indexes[i])] < node->GCX){ //West wing
                    NWi[NWc]=indexes[i];
                    NWc++;
                } else {
                    NEi[NEc]=indexes[i];
                    NEc++;
                }
            }
        }
        if(NEc>0){ //There are particles in the NE node
            node->children[0]= (Node*)malloc(sizeof *node->children[0]);
            node->children[0]->TRX=node->TRX;
            node->children[0]->TRY=node->TRY;
            node->children[0]->LLX=node->GCX;
            node->children[0]->LLY=node->GCY;
            node->children[0]->GCX=(node->GCX+node->TRX)/2;
            node->children[0]->GCY=(node->GCY+node->TRY)/2;
			//Recursively call the build Tree function in the children nodes
            buildTree(node->children[0],shrdBuff,NEi,NEc);
        } else {
            node->children[0]=NULL;
        }
        if(NWc>0){
            node->children[1]= (Node*)malloc(sizeof *node->children[1]);
            node->children[1]->TRX=node->GCX;
            node->children[1]->TRY=node->TRY;
            node->children[1]->LLX=node->LLX;
            node->children[1]->LLY=node->GCY;
            node->children[1]->GCX=(node->LLX+node->GCX)/2;
            node->children[1]->GCY=(node->GCY+node->TRY)/2;
            buildTree(node->children[1],shrdBuff,NWi,NWc);
        } else {
            node->children[1]=NULL;
        }
        if(SWc>0){
			node->children[2] = (Node*)malloc(sizeof *node->children[2]);
            node->children[2]->TRX=node->GCX;
            node->children[2]->TRY=node->GCY;
            node->children[2]->LLX=node->LLX;
            node->children[2]->LLY=node->LLY;
            node->children[2]->GCX=(node->LLX+node->GCX)/2;
            node->children[2]->GCY=(node->LLY+node->GCY)/2;
            buildTree(node->children[2],shrdBuff,SWi,SWc);
        } else {
            node->children[2]=NULL;
        }
        if(SEc>0){
            node->children[3]= (Node*)malloc(sizeof *node->children[3]);
            node->children[3]->TRX=node->TRX;
            node->children[3]->TRY=node->GCY;
            node->children[3]->LLX=node->GCX;
            node->children[3]->LLY=node->LLY;
            node->children[3]->GCX=(node->GCX+node->TRX)/2;
            node->children[3]->GCY=(node->LLY+node->GCY)/2;
            buildTree(node->children[3],shrdBuff,SEi,SEc);
        } else {
            node->children[3]=NULL;
        }
        node->mass=0;
        node->CMX=0;
        node->CMY=0;
		//Now that every children node has been created, we calculate the center of Mass
		//And total mass of the node
        for(i=0;i<4;i++){
         	   if(node->children[i]!=NULL){
                node->mass+=node->children[i]->mass;
                node->CMX+=node->children[i]->CMX*node->children[i]->mass;
                node->CMY+=node->children[i]->CMY*node->children[i]->mass;
            }
        }
        node->CMX=node->CMX/node->mass;
        node->CMY=node->CMY/node->mass;
    }
}

void calculateForce(struct Node *tree, double *shrdBuff, double *localBuff, int index){
    double distance = sqrt((tree->CMX-shrdBuff[PX(index)])*(tree->CMX-shrdBuff[PX(index)])+
                           (tree->CMY-shrdBuff[PY(index)])*(tree->CMY-shrdBuff[PY(index)]));
    if(distance>0){ //To avoid self interactions
        if(distance>rcutoff || tree->external){
            double f;
            if(distance<rlimit){
                f=G*tree->mass/(rlimit*rlimit*distance);
            } else {
                f=G*tree->mass/(distance*distance*distance);
            }
			//Projection of the acceleration in the x direction
            localBuff[AX(index)]+=f*(tree->CMX-shrdBuff[PX(index)]);
			//And in the y direction
            localBuff[AY(index)]+=f*(tree->CMY-shrdBuff[PY(index)]);
        } else { //The node is too close to the particle, recursively check each node
            int i;
            for(i=0;i<4;i++){
                if(tree->children[i]!=NULL){
                    calculateForce(tree->children[i],shrdBuff,localBuff,index);
                }
            }
        }
    }
}

void moveParticle(double *shrdBuff, double *localBuff, int index){
	//Really simple euler method to advance in time
    double oldX=shrdBuff[PX(index)];
    double oldY=shrdBuff[PY(index)];
    shrdBuff[PX(index)]+=localBuff[VX(index)]*dt+localBuff[AX(index)]*dt*dt*0.5;
    shrdBuff[PY(index)]+=localBuff[VY(index)]*dt+localBuff[AY(index)]*dt*dt*0.5;
    localBuff[VX(index)]=(shrdBuff[PX(index)]-oldX)/dt;
    localBuff[VY(index)]=(shrdBuff[PY(index)]-oldY)/dt;
}

//OpenGL Function
void drawParticle(double *shrdBuff, double *radius, int index){
    glBegin(GL_TRIANGLE_FAN);
    int k;
    glVertex2f(shrdBuff[PX(index)],shrdBuff[PY(index)]);
    for(k=0;k<20;k++){
        float angle=(float) (k)/19*2*3.141592;
        glVertex2f(shrdBuff[PX(index)]+radius[index]*cos(angle),shrdBuff[PY(index)]+radius[index]*sin(angle));
    }
    glEnd();
}

//OpenGL Function
void drawBarnesHutDivisions(struct Node *rootNode){
    if(!rootNode->external){
        glBegin(GL_LINES);
        glVertex2f(rootNode->GCX,rootNode->LLY);
        glVertex2f(rootNode->GCX,rootNode->TRY);
        glVertex2f(rootNode->LLX,rootNode->GCY);
        glVertex2f(rootNode->TRX,rootNode->GCY);
        glEnd();
        int i;
        for(i=0;i<4;i++){
            if(rootNode->children[i]!=NULL){
                drawBarnesHutDivisions(rootNode->children[i]);
            }
        }
    }
}

int main(int argc, char *argv[]){
	

	MPI_Init(&argc,&argv);
	int rank;
	int worldSize;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&worldSize);

	//nOriginal will hold the initial value of particles (some might exit the calculation area)
	int nOriginal=500;
	//nShared are the real number of particles that are being calculated
    int nShared=500;
	//Total steps for the calculation (ignored in live view)
	int steps=100;

	if(argc>1){
		nShared=atoi(argv[1]);
		if(argc>2){
			steps=atoi(argv[2]);
		}
	}

    int nLocal=ceil(((float) nShared)/worldSize);
	int nLocalMax=nLocal;
	nOriginal=nShared;
	//If the number of particles is not evenly divisible by the number of procs
	if(rank==(worldSize-1)){
		nLocal=nShared-(worldSize-1)*nLocal;
	}
	int nLocalOriginal=nLocal;

	//Buffer to hold position in x and y, and mass.
	//First entrance holds the number of particles per core
    double *sharedBuff = (double *) malloc(sizeof(double)*(3*nShared+1));
	//Buffer to receive data from MPI
    double *receiveBuff = (double *) malloc(sizeof(double)*(3*nLocalMax+1));
    //Buffer to hold velocity and acceleration in x and y directions
	double *localBuff = (double *) malloc(sizeof(double)*(4*nLocal));
	//For OpenGL, radius of the particles (constant density ~ sqrt(Mass))
    double *radius = (double *) malloc(sizeof(double)*(nShared));

	//Rand need to work different per proc
    srand(time(NULL)+rank);
    int i;
    for(i=0;i<nLocal;i++){
        sharedBuff[PX(i)]=(float) (i+rank*nLocalMax)/(nShared-1)*0.8+0.1;
        //sharedBuff[PY(i)]=(float) (rand()%4096)/4095*0.8+0.1;
        sharedBuff[PY(i)]=(float) (1+sin(i+rank*nLocalMax))/2*0.8+0.1;
        //sharedBuff[MASS(i)]=(double) (rand()%2048)/2047*2+1;
        sharedBuff[MASS(i)]=(double) (1+sin(rank*nLocalMax))/2*2+1;

        localBuff[VX(i)]=0;
        localBuff[VY(i)]=0;
        localBuff[AX(i)]=0;
        localBuff[AY(i)]=0;
    }
	sharedBuff[0]=nLocal;

	int generalCounter=nLocal;
	for(i=0;i<worldSize;i++){
		if(i==rank){
			MPI_Bcast(sharedBuff,3*nLocalMax+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
		} else {		
			MPI_Bcast(receiveBuff,3*nLocalMax+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
    		int k;
			for(k=0;k<receiveBuff[0];k++){
				sharedBuff[PX(generalCounter)]=receiveBuff[PX(k)];
				sharedBuff[PY(generalCounter)]=receiveBuff[PY(k)];
				sharedBuff[MASS(generalCounter)]=receiveBuff[MASS(k)];
        		generalCounter++;
			}
		}
	}

    struct Node* tree = (Node*)malloc(sizeof *tree);
    tree->LLX=0;
    tree->LLY=0;
    tree->TRX=1;
    tree->TRY=1;
    tree->GCX=0.5;
    tree->GCY=0.5;

    int *indexes = (int*) malloc(sizeof(int)*nShared);
    for(i=0;i<nShared;i++){
        indexes[i]=i;
    }
	
	if(rank==0){
		if(argc>=3){

			for(i=0;i<nShared;i++){
        		radius[i]=sqrt(sharedBuff[MASS(i)])*0.0025;
			}

	    	if(!glfwInit()){
    	    	//printf("Failed to start GLFW\n");
				cout << "Failed to start GLFW\n";
        		return -1;
    		}
    		GLFWwindow *window = glfwCreateWindow(2000,2000,"Simulation",NULL,NULL);
    		if(!window){
        		//printf("Failed to open window\n");
				cout << "Failed to open window\n";
        		return -1;
    		}
    		glfwMakeContextCurrent(window);
    		glfwSwapInterval(1);

    		glMatrixMode(GL_PROJECTION);
    		glLoadIdentity();
    		glOrtho(0,1,0,1,0,1);
    		glMatrixMode(GL_MODELVIEW);

    		while(!glfwWindowShouldClose(window)){
        		glClear(GL_COLOR_BUFFER_BIT);

				double t=glfwGetTime();
        		buildTree(tree,sharedBuff,indexes,nShared);
				
				nShared=nOriginal-(nLocalOriginal-nLocal);	
        		for(i=0;i<nLocal;i++){
            		localBuff[AX(indexes[i])]=0;
            		localBuff[AY(indexes[i])]=0;
            		int s;
            		for(s=0;s<4;s++){
                		if(tree->children[s]!=NULL)
                		calculateForce(tree->children[s],sharedBuff,localBuff,indexes[i]);
            		}
            		moveParticle(sharedBuff,localBuff,indexes[i]);
            		if(sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1){
                		int r;
                		nLocal--;
                		nShared--;
						sharedBuff[MASS(indexes[i])]=-1;
                		for(r=i;r<nLocal;r++){
                    		indexes[r]=indexes[r+1];
                		}
              		  	i--;
            		}
        		}

				
				int generalCounter=nLocalOriginal;
				int offset=nLocalOriginal-nLocal;
				for(i=0;i<worldSize;i++){
					if(i==rank){
						MPI_Bcast(sharedBuff,3*nLocalMax+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
					} else {		
						MPI_Bcast(receiveBuff,3*nLocalMax+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
    					int k;
						for(k=0;k<receiveBuff[0];k++){
							if(receiveBuff[MASS(k)]<0){
								offset++;
								nShared--;	
							} else {
								sharedBuff[PX(generalCounter)]=receiveBuff[PX(k)];
								sharedBuff[PY(generalCounter)]=receiveBuff[PY(k)];
								sharedBuff[MASS(generalCounter)]=receiveBuff[MASS(k)];
								indexes[generalCounter-offset]=generalCounter;
							}
        					generalCounter++;
						}
					}
				}
				

        		drawBarnesHutDivisions(tree);
        		int k;
        		for(k=0;k<nShared;k++){
            		drawParticle(sharedBuff,radius,indexes[k]);
        		}

				t=glfwGetTime()-t;
				if(t<0.013){
					Sleep(1000*1000*(0.013-t));
				}

        		glfwSwapBuffers(window);
        		glfwPollEvents();
    		}
    		glfwTerminate();
		} else {
			system("mkdir res");
			int count=0;
    		while(count<steps){
        		buildTree(tree,sharedBuff,indexes,nShared);
				
				nShared=nOriginal-(nLocalOriginal-nLocal);	
        		for(i=0;i<nLocal;i++){
            		localBuff[AX(indexes[i])]=0;
            		localBuff[AY(indexes[i])]=0;
            		int s;
            		for(s=0;s<4;s++){
                		if(tree->children[s]!=NULL)
                		calculateForce(tree->children[s],sharedBuff,localBuff,indexes[i]);
            		}
            		moveParticle(sharedBuff,localBuff,indexes[i]);
            		if(sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1){
                		int r;
                		nLocal--;
                		nShared--;
						sharedBuff[MASS(indexes[i])]=-1;
                		for(r=i;r<nLocal;r++){
                    		indexes[r]=indexes[r+1];
                		}
              		  	i--;
            		}
        		}

				
				int generalCounter=nLocalOriginal;
				int offset=nLocalOriginal-nLocal;
				for(i=0;i<worldSize;i++){
					if(i==rank){
						MPI_Bcast(sharedBuff,3*nLocalMax+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
					} else {		
						MPI_Bcast(receiveBuff,3*nLocalMax+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
    					int k;
						for(k=0;k<receiveBuff[0];k++){
							if(receiveBuff[MASS(k)]<0){
								offset++;
								nShared--;	
							} else {
								sharedBuff[PX(generalCounter)]=receiveBuff[PX(k)];
								sharedBuff[PY(generalCounter)]=receiveBuff[PY(k)];
								sharedBuff[MASS(generalCounter)]=receiveBuff[MASS(k)];
								indexes[generalCounter-offset]=generalCounter;
							}
        					generalCounter++;
						}
					}
				}
				
				char filename[]={'r','e','s','/','r','0','0','0','0'};
				filename[8]=count%10+48;
				filename[7]=(count/10)%10+48;
				filename[6]=(count/100)%10+48;
				filename[5]=(count/1000)%10+48;
				ofstream res;
				res.open(string(filename));
				//FILE *res = fopen(filename,"w");
				for(i=0;i<nShared;i++){
					//fprintf(res,"%d\t%e\t%e\n",indexes[i],sharedBuff[PX(indexes[i])],sharedBuff[PY(indexes[i])]);
					res << indexes[i] << sharedBuff[PX(indexes[i])] << sharedBuff[PY(indexes[i])];
					res << "\n";
				}
				//fclose(res);
				res.close();
				count++;
			}	
		}
	} else {
		int count=0;
		int dec=1;
		if(argc>3){
			dec=0;
		}
    	while(count<steps){
        	buildTree(tree,sharedBuff,indexes,nShared);
			
			nShared=nOriginal-(nLocalOriginal-nLocal);	
        	for(i=0;i<nLocal;i++){
            	localBuff[AX(indexes[i])]=0;
            	localBuff[AY(indexes[i])]=0;
            	int s;
            	for(s=0;s<4;s++){
                	if(tree->children[s]!=NULL)
                	calculateForce(tree->children[s],sharedBuff,localBuff,indexes[i]);
            	}
            	moveParticle(sharedBuff,localBuff,indexes[i]);
            	if(sharedBuff[PX(indexes[i])]<=0 || sharedBuff[PX(indexes[i])]>=1 || sharedBuff[PY(indexes[i])] <=0 || sharedBuff[PY(indexes[i])] >= 1){
                	int r;
                	nLocal--;
                	nShared--;
					sharedBuff[MASS(indexes[i])]=-1;
                	for(r=i;r<nLocal;r++){
                    	indexes[r]=indexes[r+1];
                	}
       		  	i--;
				}
       		}
			int generalCounter=nLocalOriginal;
			int offset=nLocalOriginal-nLocal;
			for(i=0;i<worldSize;i++){
				if(i==rank){
					MPI_Bcast(sharedBuff,3*nLocalMax+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
				} else {		
					MPI_Bcast(receiveBuff,3*nLocalMax+1,MPI_DOUBLE,i,MPI_COMM_WORLD);
    				int k;
					for(k=0;k<receiveBuff[0];k++){
						if(receiveBuff[MASS(k)]<0){
							offset++;
							nShared--;	
						} else {
							sharedBuff[PX(generalCounter)]=receiveBuff[PX(k)];
							sharedBuff[PY(generalCounter)]=receiveBuff[PY(k)];
							sharedBuff[MASS(generalCounter)]=receiveBuff[MASS(k)];
							indexes[generalCounter-offset]=generalCounter;
						}
        				generalCounter++;
					}
				}
			}
			count+=dec;
        }
	}
	
	MPI_Finalize();
    return 0;
}