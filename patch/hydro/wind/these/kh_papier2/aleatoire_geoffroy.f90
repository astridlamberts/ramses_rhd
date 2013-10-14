double randm(void) {
const int a = 16807;
const int m = 2147483647;
static int in0 = 13763;
int q;


// When using mpi, this allows us to have different series in each process...

    
#ifdef PARALLEL
if(in0 == 13763) in0 += 2543 * prank;
#endif

    
/* find random number  */
q= (int) fmod((double) a * in0, m);
in0=q;


return((double)q/(double)m);
}
