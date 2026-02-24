#include "PPICLF_STD.h"
#:include 'PPICLF_PARTMACROS.fypp'

module ppiclf_m_particledata
    use ppiclf_m_types
    implicit none

    ! from [-PPICLF_LPART_GP to 0) is for ghost particles recieved from other ranks
    ! [0,0] is reserved for temporary particles as needed (for periodicity I think)
    ! (0, PPICLF_LPART] is for real particles owned by this rank
    @:DECLAREPARTVAR(PPICLF_t_particle, ppiclf_parts, (PPICLF_LPART))
    @:DECLAREPARTVAR(PPICLF_t_ghostParticle, ppiclf_gparts, (PPICLF_LPART_GP))

    contains

    subroutine CopyRealToGhost(@{LISTCOMPONENTS(PPICLF_t_particle, particle)}@, @{LISTCOMPONENTS(PPICLF_t_ghostParticle, ghost)}@)
        @:DECLAREPARTVAR(PPICLF_t_particle, particle)
        @:DECLAREPARTVAR(PPICLF_t_ghostParticle, ghost)
        

    end subroutine CopyRealToGhost
end module ppiclf_m_particledata