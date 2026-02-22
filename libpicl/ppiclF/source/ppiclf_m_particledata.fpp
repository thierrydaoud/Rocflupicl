#include "PPICLF_STD.h"
module ppiclf_m_particledata
    use ppiclf_m_types, only : ppiclf_t_particlepos, ppiclf_t_particle
    implicit none

    ! from [-PPICLF_LPART_GP to 0) is for ghost particles recieved from other ranks
    ! [0,0] is reserved for temporary particles as needed (for periodicity I think)
    ! (0, PPICLF_LPART] is for real particles owned by this rank
    type(PPICLF_t_particlepos) :: ppiclf_partpos(-PPICLF_LPART_GP:PPICLF_LPART)
    type(PPICLF_t_particle)     :: ppiclf_parts(-PPICLF_LPART_GP:PPICLF_LPART)
end module ppiclf_m_particledata