program timcom_omip
  implicit none

#ifdef USE_GPU
  call timcom_main_gpu
#else
  call timcom_main
#endif

end program timcom_omip
