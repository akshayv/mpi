#include <stdio.h>
#include"mpi.h"
#include<stdlib.h>
#include<math.h>
#include<time.h>

int SIZE = 32;
int NUMBER_OF_STICKS = 16;

int DIAMETER_THRESHOLD = 6;
float AVERAGE_DIST_THRESHOLD = 3.2984;

int MAX_NUMBER_ITER = 100;

int world_rank;

typedef struct Node {
  int id;
  int leftNeighbourId;
  int middleNeighbourId;
  int rightNeighbourId;
} Node;

typedef struct DiameterDetails {
  int diameter;
  float sumDistance;
} DiameterDetails;

Node* alloc_node_array(int number);
int* alloc_int_array(int number);
void generateConnectingSticks(Node* nodes, int numberOfSticks);
void deleteNodeAt(int* nodes, int index, int size);
DiameterDetails getDiameterDetails(Node* nodes, int startNode, int endNode);
int getMinDistance(Node node, Node secondNode, int currentDistance, Node* nodes);
Node* findById(Node* nodes, int id, int size);
void print(Node* nodes, int size);

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  //int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size < 2) {
           fprintf(stderr, "World size must be greater than 1");
           MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  int modulo = SIZE % world_size;
  if(modulo != 0) {
      SIZE = SIZE + (SIZE - ((SIZE / world_size) * world_size));
  }

  int nodesPerCore = SIZE/world_size;
  
  Node* nodes;
  nodes = alloc_node_array(SIZE);
  int i;
  int number_iterations = 0;
  double total_time = 0;
  double t1, t2; 
  t1 = MPI_Wtime();
  
  if(world_rank == 0) {
    for(i = 0; i < SIZE; i++) {
      nodes[i].id = i;
      nodes[i].leftNeighbourId = getMod(i + 1, SIZE);
      nodes[i].rightNeighbourId = getMod(i - 1, SIZE);
    }
  }

  while(1) {
    number_iterations++;
    if(world_rank == 0) {
      generateConnectingSticks(nodes, NUMBER_OF_STICKS); 
    } 
    
    MPI_Bcast(nodes, SIZE * sizeof(Node), MPI_BYTE, 0, MPI_COMM_WORLD);
    //t1 = MPI_Wtime();
    DiameterDetails diameterDetails;
    int startNode = world_rank * nodesPerCore;
    int endNode = (world_rank + 1) * nodesPerCore;
    diameterDetails = getDiameterDetails(nodes, startNode, endNode);
    //t2 = MPI_Wtime();
    //total_time += (t2-t1)/100;
    int globalMax = 0;
    float globalSum;
    MPI_Reduce(&(diameterDetails.diameter), &(globalMax), 1,
             MPI_INT, MPI_MAX, 0,
             MPI_COMM_WORLD);
  
    MPI_Reduce(&(diameterDetails.sumDistance), &(globalSum), 1,
             MPI_FLOAT, MPI_SUM, 0,
             MPI_COMM_WORLD);
    
    int status = 0;
    if( world_rank == 0 ) {
      float globalAvg = globalSum / ((SIZE) * (SIZE -1));
      printf("Diameter: %d, AverageDistance: %f \n", globalMax, globalAvg);
      fflush(stdout);
      if(globalMax <= DIAMETER_THRESHOLD 
        && globalAvg <= AVERAGE_DIST_THRESHOLD) {
	DIAMETER_THRESHOLD = globalMax;
	AVERAGE_DIST_THRESHOLD = globalAvg;
        print(nodes, SIZE);
      } 
      if(number_iterations >= MAX_NUMBER_ITER) {
        status = 1;
      }
    }
    
    
    MPI_Bcast(&status, sizeof(int), MPI_INT, 0, MPI_COMM_WORLD);
    if(status == 1) {
      break;
    }
    
  }
  
  t2 = MPI_Wtime();
  total_time = t2-t1;
  printf( "Elapsed time is %f\n", total_time ); 
  printf( "Average time is %f\n", (total_time)/number_iterations ); 
      
  MPI_Finalize();
}

int getMod(int a, int b){ return (a%b+b)%b; }

Node* alloc_node_array(int number) {
  Node *data = (Node *)calloc(number, sizeof(Node));
  return data;
}

int* alloc_int_array(int number) {
  int *data = (int *)calloc(number, sizeof(int));
  return data;
}

void print(Node* nodes, int size) {
  int i;
  for(i = 0; i < size; i++) {
    printf("%d -> %d, %d, %d", nodes[i].id, nodes[i].leftNeighbourId,
	   nodes[i].middleNeighbourId, nodes[i].rightNeighbourId);
    printf("\n");
  }
  printf("\n");
}

void generateConnectingSticks(Node* nodes, int numberOfSticks) {
  int* availableNodes;
  availableNodes = alloc_int_array(SIZE);
  int i; 
  srand(time(NULL) + rand());
  for(i = 0; i < SIZE; i++) {
    availableNodes[i] = nodes[i].id;
  }
  int size = SIZE;
  for(i = 0; i < numberOfSticks; i++) {
    Node* node = findById(nodes, availableNodes[0], SIZE);
    deleteNodeAt(availableNodes, 0, size--);
    int r = rand() % (size);
    Node* connectedNode = findById(nodes, availableNodes[r], SIZE);
    node->middleNeighbourId = connectedNode->id;
    connectedNode->middleNeighbourId = node->id;
    deleteNodeAt(availableNodes, r, size--);
  }
}

void deleteNodeAt(int* nodes, int index, int size) {
  int i;
  for(i = index; i < size - 1; i++) {
    nodes[i] = nodes[i + 1];
  }
}

Node* findById(Node* nodes, int id, int size) {
  int i;
  for(i = 0; i < size; i++) {
    if(id == nodes[i].id) return &(nodes[i]);  
  }
}

DiameterDetails getDiameterDetails(Node* nodes, int startNode, int endNode) {
  int maxDiameter = -1;
  float sumDistance = 0;
  int currMin;
  int i, j;   
  for(i = startNode; i < endNode; i++) {
    for(j = 0; j < SIZE; j++) {
      if(i == j) continue;
      currMin = getMinDistance(nodes[i], nodes[j], 0, nodes);
      sumDistance += currMin;
      if (currMin > maxDiameter) {
	maxDiameter = currMin;
      }
    }
  }
  
  DiameterDetails details = {.diameter = maxDiameter, .sumDistance = sumDistance};
  return details;
}

int getMinDistance(Node node, Node secondNode, int currentDistance, Node* nodes) {
  if(secondNode.id == node.id)  return currentDistance;
  int minDistance = -1;
  int currMin;
  if (currentDistance > (SIZE/4)) return minDistance;
  currMin = getMinDistance(nodes[node.leftNeighbourId], secondNode, currentDistance + 1, nodes);
  if (currMin != -1)
    minDistance = currMin;
  currMin = getMinDistance(nodes[node.middleNeighbourId], secondNode, currentDistance + 1, nodes);
  if ((currMin < minDistance || minDistance == -1) && currMin != -1)
    minDistance = currMin;
  currMin = getMinDistance(nodes[node.rightNeighbourId], secondNode, currentDistance + 1, nodes);
  if ((currMin < minDistance || minDistance == -1) && currMin != -1)
    minDistance = currMin;
  return minDistance;
}


