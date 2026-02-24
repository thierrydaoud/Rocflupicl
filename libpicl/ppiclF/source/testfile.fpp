#:include 'PPICLF_PARTMACROS.fypp'



@:GENERATESTRUCTS()



@:DECLAREPARTVAR(PPICLF_t_particle, ppiclf_parts, (PPICLF_LPART))


@{USEPARTICLE(ppiclf_parts(i)%y%pos%x)}@

@{USEPARTICLE(ppiclf_parts(i)%yodtc%vel%z)}@

@{USEPARTICLE(ppiclf_parts(i, 1:10)%y%pos%x)}@
