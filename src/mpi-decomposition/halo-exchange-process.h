#ifndef ONZA_FDTD_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
#define ONZA_FDTD_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
///
/// @file   halo-exchange-process.h
/// @author Ladutenko Konstantin <kostyfisik at gmail (.) com>
/// @date   Thu Apr 26 19:27:05 2012
/// 
/// @brief  Exchange borders of computational domain 

/// @brief Class for MPI process.
/// Contains computational domain borders data, methods to exchange
/// and update borders.
class HaloExchangeProcess {
 public:
  /// @brief Accesor
  ///
  /// @return associated with object rank of the process
  int process_rank() const { return process_rank_ ; }
  int init();
 private:
  int process_rank_;
};
#endif // ONZA_FDTD_MPI_DECOMPOSITION_HALO_EXCHANGE_PROCESS_H_
