#include "libinclude.h"
#include "const.h"
#include "fundec.h"

vector<double> realnreci0(){
    // This term will contain the real energy part as energy[0] and reciprocal energy k=0 as energy[1]
    vector<double> energy(0,2);
    omp_set_num_threads(thread::hardware_concurrency());
    for (int i = 0; i < natoms; i++){
        for (int j = 0; j < i; j++){
            
        }
        
    }
    
    // real energy calculation


    return out;
}