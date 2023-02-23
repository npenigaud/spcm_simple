module cublas_mod
!
! Define the INTERFACE to the NVIDIA C code cublasSgemm and cublasDgemm
!
interface cuda_gemm
!
! void cublasSgemm (char transa, char transb, int m, int n,
! int k, float alpha, const float *A, int lda,
! const float *B, int ldb, float beta, float *C, int ldc)
!
subroutine cuda_sgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='cublasSgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_float),value :: alpha,beta
real(c_float), device, dimension(lda,*) :: A
real(c_float), device, dimension(ldb,*) :: B
real(c_float), device, dimension(ldc,*) :: C
end subroutine cuda_sgemm

!
! void cublasDgemm (char transa, char transb, int m, int n,
! int k, double alpha, const double *A, int lda,
! const double *B, int ldb, double beta, double *C, int ldc)
!
subroutine cuda_dgemm(cta, ctb, m, n, k,&
alpha, A, lda, B, ldb, beta, c, ldc) bind(C,name='cublasDgemm')
use iso_c_binding
character(1,c_char),value :: cta, ctb
integer(c_int),value :: m,n,k,lda,ldb,ldc
real(c_double),value :: alpha,beta
real(c_double), device, dimension(lda,*) :: A
real(c_double), device, dimension(ldb,*) :: B
real(c_double), device, dimension(ldc,*) :: C
end subroutine cuda_dgemm

end interface

end module cublas_mod

