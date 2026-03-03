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

#:for n, real_ref, real_off, ghost_ref, ghost_off in fyppmacros.Loop_All_Reals("particle", "ghost", sameArrays=True)
        ${ghost_ref}$(${ghost_off + 1}$: ${ghost_off + n}$)  = ${real_ref}$(${real_off + 1}$: ${real_off + n}$)
#:endfor

        @{USEPARTICLE(ghost%iprop)}@ = @{USEPARTICLE(particle%iprop)}@
    end subroutine CopyRealToGhost
end module ppiclf_m_particledata