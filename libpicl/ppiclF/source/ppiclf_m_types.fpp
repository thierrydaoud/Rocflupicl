#:include 'PPICLF_PARTMACROS.fypp'

module ppiclf_m_types

    ! make the particle id a structure so that it could hold more info in the future if needed
    ! ie, id could be 3 parts, a node id, a timestep # and a sequential number, to ensure that every particle has a unique ID
    ! since the ID follows the particle as it moves between ranks, it would allow for the particle to be tracked back to when(timestep) and where(node) it was added to the simulation
    ! for now though its just a single number for compatability
    type :: ppiclf_part_id
        integer global_id
    end type ppiclf_part_id

    @:GENERATESTRUCTS()


end module ppiclf_m_types




