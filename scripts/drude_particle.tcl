#To call first:

#inter $bondId_drude drude $temp_core $gamma_core $temp_drude $gamma_drude $k_drude $mass_red_drude $r_cut

proc add_drude_to_core { bondId_drude id_core id_drude type_drude polarization initialSpread } {

	#PARAMETERS:
	#--------------------------------------------------------------------------------
	#id_core:	particle ID of existing core particle 
	#id_drude:	free particle ID for drude particle
	#type_drude:	particle type for drude particle
	#sigma_core:	core particle 'size' to recalculate polarization
        #polarization:  in [length]^-3 gets unitless with particle volume (*sigma_core^-3)

	#SCRIPT-DESCRIPTION:
	#--------------------------------------------------------------------------------
	#		-Adds drude particle to existing particle 'id_core'
	#		-Disables thermostat for core and drude via LPP and temp=gamma=0
	#		-Adds drude-bond between core and drude: harmonic bond + langevin on relative coordinates
	#		-Adds subt_elec-bond between core and drude: Subtracts electrostatic interactino 
        #		-Splits mass between core and drude
        #		-Splits charge between core and drude (see formula for q_drude)

	#REQUIREMENTS:
	#--------------------------------------------------------------------------------
	#		Features: LANGEVIN_PER_PARTICLE, ELECTROSTATICS, MASS
	#		Existing charged core particle
	#		Initialization call 'init_drude'

	#global bondId_drude_4s3jkhdf82hka
	#global bondId_subt_elec_hdfkj4hrkjfy
	#global k_drude_skjfhes34hkjss7

	#if {$bondId_drude_4s3jkhdf82hka == ""} {
	#	return "ERROR: Call init_drude first to set bond IDs."
	#} elseif {[part $id_core print q] == 0} {
	#	return "ERROR: Core particle has to be charged."
	#}

	set warnings ""

	if {$polarization <= 0} {
		return "ERROR: Polarization must be a positive number."
        } elseif { [part $id_core print] == "na" } {
		return "ERROR: Can't find core particle with id $id_core."
        } elseif { [inter $bondId_drude] == "unknown bonded interaction number $bondId_drude" } {
		return "ERROR: Can't find drude bond with bond id $bondId_drude. Create a drude-bond interaction type first."
        } elseif { [part $id_drude] != "na" } {
		set warnings "WARNING: Particle with id $id_drude already exists."
        }

	set gamma_core [lindex [inter $bondId_drude] 4]
	set gamma_drude [lindex [inter $bondId_drude] 6]
	set k_drude [lindex [inter $bondId_drude] 7]
	set mass_core [part $id_core print mass]
	set mass_drude [expr $mass_core * [lindex [inter $bondId_drude] 8]]
	set q_core [part $id_core print q]

	if {$q_core > 0} {
		set sign_q 1.0
	} else {
 		set sign_q -1.0
        }
#	set q_drude [expr $sign_q * pow($k_drude*$polarization/pow($sigma_core,3.0), 1./2.)]
	set q_drude [expr $sign_q * pow($k_drude*$polarization, 1./2.)]

	
	part $id_core mass [expr $mass_core-$mass_drude] q [expr $q_core - $q_drude] temp 0 gamma 0
	set drude_px [expr [lindex [part $id_core print pos] 0] + 2.0*([t_random]-0.5)*$initialSpread]
	set drude_py [expr [lindex [part $id_core print pos] 1] + 2.0*([t_random]-0.5)*$initialSpread]
	set drude_pz [expr [lindex [part $id_core print pos] 2] + 2.0*([t_random]-0.5)*$initialSpread]
	part $id_drude pos $drude_px $drude_py $drude_pz v 0 0 0 q $q_drude type $type_drude mass $mass_drude temp 0 gamma 0
	part $id_core bond $bondId_drude $id_drude
	#part $id_core bond $bondId_subt_elec_hdfkj4hrkjfy $id_drude

	return "Drude particle created with parameters:\n\
		id: $id_drude\n\
		charge: $q_drude\n\
		mass: $mass_drude\n\
		Core charge changed to: [expr $q_core - $q_drude]\n\
		Core mass changed to: [expr $mass_core-$mass_drude]\n\
		Relaxation time drude thermostat: [expr $mass_drude/$gamma_drude]\n\
		Period spring: [expr 2.0*[PI]*sqrt($mass_drude/$gamma_drude)]\n\
		Relaxation time core thermostat: [expr ($mass_core-$mass_drude)/$gamma_core]\n\
		$warnings"
}
