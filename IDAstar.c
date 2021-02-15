// -*- C++ -*-  (hint: gcc -Wall  -o $(PGM) $(PGM).c -I.)
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
/**@file  IDAstar.c
 * @brief solve 4x4 sliding puzzle with IDA* algorithm
 *
 *  https://algorithmsinsight.wordpress.com/graph-theory-2/ida-star-algorithm-in-general/
 * see https://en.wikipedia.org/wiki/Iterative_deepening_A*
 * 2021-feb-15
 * by RMM
 */ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef unsigned char u8t;
typedef unsigned short u16t;
enum { NR=4, NC=4, NCELLS = NR*NC };
enum { RIGHT, LEFT, UP, DOWN, NDIRS };
enum { OK = 1<<8, XX = 1<<9, FOUND = 1<<10 };
enum { MAX_INT=0x7E, MAX_NODES=(16*65536)*400};
enum { BIT_HDR=1<<0, BIT_GRID=1<<1, BIT_OTHER=1<<2 };
typedef struct { u16t dn; u16t hn; }HSORT_T;

typedef struct { unsigned id; unsigned src; u8t data[NCELLS];
   u8t h; u8t g; u8t udlr;}NODE_T;

NODE_T goal={0,0,{1,2,3,4, 5,6,7,8, 9,10,11,12, 13,14,15,0},0,0,0};
NODE_T work; // copy of puzzle with changes

NODE_T G19={0,0, //g=19,n=26
	    {9,5,1,4, 10,6,2,8,  13,14,3,11, 0,15,7,12},0,0,0};

NODE_T G28={0,0, //g=28; n=1261
	    {9,5,1,4, 13,10,6,8, 15,3,2,11, 14,0,7,12},0,0,0};

NODE_T G34={0,0, //g=34; n=248,055; seconds=0.6
	    {13,9,5,4, 15,6,1,8, 0,10,2,11, 14,3,7,12},0,0,0};

NODE_T G52={0,0, //(g=52; n=245M;) n=34,296,567; minutes=1.3
	    {15,13,9,5, 14,6,1,4,  10,12,0,8, 3,7,11,2},0,0,0};

NODE_T G53={0,0, //(g=52+; fail;) g=53; n=116,346,763; minutes=4.5
	    {15,13,9,5, 14,6,0,4, 10,12,1,8, 3,7,11,2},0,0,0};

NODE_T G99={0,0, // formidable; TBD; should be at G=52
	    {15,14,1,6, 9,11,4,12, 0,10,7,3, 13,8,5,2},0,0,0};

struct {
   unsigned nodes;
   int stack_idx;
   int max_stack;
   u16t gfound;
   unsigned root_visits;
   unsigned verbose;
}my;

u16t  IDA_star(NODE_T *pNode);
u16t  make_node(NODE_T *pNode, NODE_T *pNew, u8t udlr );
u16t  search(NODE_T *pNode, u16t bound);
u16t  taxi_dist( NODE_T *pNode, u8t g );
u16t  tile_home( NODE_T *p44);
void  print_node( NODE_T *pN, const char *pMsg, short force );
u16t  goal_found(NODE_T *pNode);
char  udlr_to_char( char udlr );
void  idx_to_rc( u16t idx, u16t *row, u16t *col );
void  sort_nodes(HSORT_T *p);

//-------------------------------------------------------------------
int main( )
{
   printf( "######## IDAstar 2021-feb-15 ########\n");
   my.max_stack = 0;
   my.verbose = 0;		// minimal print node
   //my.verbose |= BIT_HDR;	// node header
   //my.verbose |= BIT_GRID;	// node grid

   memcpy(&work,&G34,sizeof(NODE_T));
   IDA_star(&work);
   return 0;
}

//-------------------------------------------------------------------
// driver for Iterative Deepining A*
u16t IDA_star(NODE_T *pN){
   my.nodes = 1;
   my.gfound = 0;
   my.root_visits = 0;
   pN->udlr = NDIRS;
   pN->g = 0;
   pN->h = taxi_dist(pN,pN->g);
   pN->id = my.nodes;
   pN->src = 0;
   const char *pr = {"Start"}; // for g++
   print_node( pN,pr,1 );
   u16t depth = pN->h;
   while(1){
      depth = search(pN,depth);
      if( depth & FOUND){
         return FOUND;  // goodbye
      }
      if( depth & 0xFF00 ){
	 printf("..error %x\n",depth);
	 return depth;
      }
      my.root_visits++;
      printf("[root visits: %u, depth %u]\n",my.root_visits,depth);
   }
   return 0;
}

//-------------------------------------------------------------------
/// search is recursive. nodes are instance variables
u16t search(NODE_T *pN, u16t bound){ 
   if(bound & 0xff00){ return bound; }
   u16t f = pN->g + pN->h;
   if( f > bound){ return f; }
   if(goal_found(pN)){
      my.gfound = pN->g;
      print_node(pN,"found",1);
      printf("nodes = %d, g=%u\n", my.nodes, my.gfound);
      return FOUND;
   }
   NODE_T news;
   // Sort successor nodes so that the lowest heuristic is visited
   // before the less promising at the same level. This reduces the
   // number of searches and finds more solutions
   HSORT_T hlist[NDIRS];
   for( short i=0; i<NDIRS; i++ ){
      u16t rv = make_node(pN,&news, i );
      hlist[i].dn = i;
      if( rv & OK ){
	 hlist[i].hn = news.h;
	 continue;
      }
      hlist[i].hn = XX;
   }
   sort_nodes(&hlist[0]);
   
   u16t temp, min = MAX_INT;
   for( short i=0; i<NDIRS; i++ ){
      if( hlist[i].hn > 0xff ) continue;
      temp = make_node(pN,&news, hlist[i].dn );
      if( temp & XX ) return XX;
      if( temp & OK ){
	 news.id = my.nodes++;
	 print_node(&news," succ",0 );
	 temp = search(&news, bound);
	 if(temp & 0xff00){  return temp;}
	 if(temp < min){ min = temp; }
      }
   }
   return min;
}

//-------------------------------------------------------------------
void sort_nodes(HSORT_T *p){
   //printf("h0=%u,h1=%u,h2=%u,h3=%u\n", p[0].hn,p[1].hn, p[2].hn,p[3].hn);
   for( short s=0; s<NDIRS-1; s++ ){
      HSORT_T tmp = p[0];
      if( p[1].hn < p[0].hn ){tmp=p[0]; p[0]=p[1]; p[1]=tmp; }
      if( p[2].hn < p[1].hn ){tmp=p[1]; p[1]=p[2]; p[2]=tmp; }
      if( p[3].hn < p[2].hn ){tmp=p[2]; p[2]=p[3]; p[3]=tmp; }
      //printf("h0=%u,h1=%u,h2=%u,h3=%u\n",p[0].hn,p[1].hn,p[2].hn,p[3].hn);
   }
}

//-------------------------------------------------------------------
/// make successor node with blank tile moved UDRL
/// return success or error
u16t make_node(NODE_T *pSrc, NODE_T *pNew, u8t udlr ){
   u16t row,col,home_idx;
   if(udlr>=NDIRS||udlr<0 ){ printf("invalid udlr %u\n",udlr); return XX; }
   if(my.nodes > MAX_NODES ){ printf("excessive nodes %u\n",my.nodes);
      return XX; }
   memcpy(pNew,pSrc,sizeof(NODE_T));
   home_idx = tile_home(pNew);
   idx_to_rc(home_idx, &row, &col );

   if( udlr == UP ){    if(row < 1)         return 0; row -= 1;  }
   if( udlr == DOWN ){  if(row >= (NR-1) )  return 0; row += 1; }
   if( udlr == LEFT){   if( col < 1 )      return 0; col -= 1; }
   if( udlr == RIGHT ){ if( col >= (NC-1) ) return 0; col += 1; }
   
   int idx2 = row * NR + col;
   if( idx2 < NCELLS ){
      u8t *p = &pNew->data[0];
      p[home_idx] = p[idx2];
      p[idx2]     = 0; // swap
      pNew->src   = pSrc->id;
      pNew->g     = pSrc->g + 1;
      pNew->h     = taxi_dist(pNew,pNew->g);
      pNew->udlr  = udlr; // latest move;
      return OK;
   }
   return 0;
}

//-------------------------------------------------------------------
/// return (HOME) position or index of blank tile
u16t tile_home(NODE_T *pN ){
   for( short i=0; i<NCELLS; i++ ){
      if( pN->data[i] == 0 ) return i;
   }
   return XX;
}

//-------------------------------------------------------------------
void print_node( NODE_T *pN, const char *pMsg, short force ){
   const int tp1 = 0;
   if( my.verbose & BIT_HDR || force || tp1){
      char ch = udlr_to_char(pN->udlr);
      printf("id:%d src:%d; h=%d, g=%u, udlr=%c, %s\n",
	     pN->id, pN->src, pN->h, pN->g, ch, pMsg);
   }
   if(my.verbose & BIT_GRID || force || tp1){
      for(u16t i=0; i<NR; i++ ){
	 for( u16t j=0; j<NC; j++ ){
	    printf("%3d",pN->data[i*NR+j]);
	 }
	 printf("\n");
      }
      printf("\n");
   }
}

//-------------------------------------------------------------------
/// return true if the new NxN equals goal
u16t goal_found(NODE_T *pN1 ) {
   for( short i=0; i<NCELLS; i++ ){
      if( pN1->data[i] != goal.data[i] ) return 0;
   }
   return 1;
}

//-------------------------------------------------------------------
/// convert UDLR index to printable char
char udlr_to_char( char udlr ){
   const char *pc = "RLUD";
   char ch = '?';
   if( udlr >= 0 && udlr < NDIRS ) ch = pc[(int)udlr];
   return ch;
}

//-------------------------------------------------------------------
/// \brief convert 1-D array index to 2-D row-column
void idx_to_rc( u16t idx, u16t *row, u16t *col ){
   *row = idx/NR; *col = abs( idx - (*row * NR));
}

//-------------------------------------------------------------------
/// \brief sum of 'manhattan taxi' distance between tile locations
u16t taxi_dist( NODE_T *pN, u8t g ){
   u8t *p44 = &pN->data[0];
   u16t tile,sum = 0;
   if(1){
      u16t r1,c1,r2,c2;
      for( u16t i=0; i<NCELLS; i++ ){
	 tile = p44[i];
	 if( tile==0 ) continue;
	 idx_to_rc(tile-1, &r1, &c1 );
	 idx_to_rc(i, &r2, &c2 );
	 sum += abs(r1-r2) + abs(c1-c2);
	 //printf("sum %d, r1=%d, c1=%d, r2=%d, c2=%d\n",
	 //sum, r1,c1, r2,c2);
      }
   }
   if(1){ // Linear Conflict by row TBD
      u8t *p44 = &pN->data[0];
      u8t a,b,c,d,L,R;
      u8t a_, b_, c_, d_;		// flags
      L=1; R=4;
      for(short row=0; row<NR; row++){
	 a=p44[0]; b=p44[1]; c=p44[2]; d=p44[3];
	 a_ = ( a>=L && a<=R );
	 b_ = ( b>=L && b<=R );
	 c_ = ( c>=L && c<=R );
	 d_ = ( d>=L && d<=R );
	 if( (a_ && b_) && a>b ) sum+=2;
	 if( (a_ && c_) && a>c ) sum+=2;
	 if( (a_ && d_) && a>d ) sum+=2;
	 if( (b_ && c_) && b>c ) sum+=2;
	 if( (b_ && d_) && b>d ) sum+=2;
	 if( (c_ && d_) && c>d ) sum+=2;
	 L+=NR; R+=NR;
      }
   }
   if(0){ // Linear Conflict by column TBD
      u8t *p44 = &pN->data[0];
      u8t a,b,c,d;
      u8t a_, b_, c_, d_;		// flags
      for(short col=0; col<NC; col++ ){
	 a = p44[col+0];
	 b = p44[col+4];
	 c = p44[col+8];
	 d = p44[col+12];
	 a_ = ( a>=1 && a<=4 );
	 b_ = ( b>=5 && b<=8 );
	 c_ = ( c>=9 && c<=12 );
	 d_ = ( d>=13 && d<=15 );
	 if( (a_ && b_) && a>b ) sum+=2;
	 if( (a_ && c_) && a>c ) sum+=2;
	 if( (a_ && d_) && a>d ) sum+=2;
	 if( (b_ && c_) && b>c ) sum+=2;
	 if( (b_ && d_) && b>d ) sum+=2;
	 if( (c_ && d_) && c>d ) sum+=2;
      }
   }

   return sum;
}

// end IDAstar.c
