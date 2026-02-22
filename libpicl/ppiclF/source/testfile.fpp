#:include 'PPICLF_PARTMACROS.fypp'



@:GENERATESTRUCTS()


type(PPICLF_t_particlepos) :: ppiclf_partpos(1000)
type(PPICLF_t_particle)     :: ppiclf_parts(1000)


@{USEPARTICLE(500, ydot, x)}@

DO i=1, 1000
#:for structArray, memberArray, n in fyppmacros.Loop_All_RealArrays()
    DO j=1, ${n}$
        ${structArray}$(i)%${memberArray}$ = 0.0
    END DO
#:endfor
#:for structArray, memberArray, n in fyppmacros.Loop_All_IntArrays()
    DO j=1, ${n}$
        ${structArray}$(i)%${memberArray}$ = 0.0
    END DO
#:endfor
END DO