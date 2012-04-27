///
/// @file   halo-exchange-process.cc
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Thu Apr 26 19:13:58 2012
/// 
/// @brief  Exchange borders of computational domain 
/// 
/// 
///
#include <cstdio>
#include <mpi.h>
#include "halo-exchange-process.h"
/// @brief Initialize process variables
/// 
///
/// @return 0 or error code.
///
int  HaloExchangeProcess::init() {  
  process_rank_ = MPI::COMM_WORLD.Get_rank( );
  //  printf("My rank is %i\n", process_rank_);
  return 0;
}
    
