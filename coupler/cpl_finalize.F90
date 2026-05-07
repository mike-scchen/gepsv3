subroutine cpl_finalize
  use m_MCTWorld, only: MCTWorld_clean => clean
  implicit none

  call MCTWorld_clean
end subroutine cpl_finalize
