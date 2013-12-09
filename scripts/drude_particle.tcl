#proc init_drude_bond { bondId_drude bondId_subt_elec temp_core gamma_core {k_drude 192.73} {temp_drude [expr $temp_core/300.0]} {gamma_drude [expr $gamma_core*15.0]} {mass_red_drude 0.1} } {

	#PARAMETERS:
	#--------------------------------------------------------------------------------
	#bondId_drude:		?
	#bondId_subt_elec:	?
	#temp_core:		?
	#gamma_core:		?
        #k_drude:  		?
	#temp_drude:		?
	#gamma_drude:		?
	#mass_red_drude:	?

#	if {[string match "*LANGEVIN_PER_PARTICLE*" [code_info]] == 0 || [string match "*ELECTROSTATICS*" [code_info]] == 0 || [string match "*MASS*" [code_info]] == 0} {
#		return "ERROR: Drude particle requires LANGEVIN_PER_PARTICLE, ELECTROSTATICS and MASS."
#        }

#	global bondId_drude_4s3jkhdf82hka
#	global bondId_subt_elec_hdfkj4hrkjfy
#	global k_drude_skjfhes34hkjss7
#	global mass_red_drude_sdhk873iassl0

#	set mass_red_drude_sdhk873iassl0 $mass_red_drude
#	set k_drude_skjfhes34hkjss7 $k_drude
#	set bondId_drude_4s3jkhdf82hka $bondId_drude
#	set bondId_subt_elec_hdfkj4hrkjfy $bondId_subt_elec

#        inter $bondId_drude drude $k_drude -1 $temp_core $gamma_core $temp_drude $gamma_drude
#        inter $bondId_subt_elec subt_elec
#}

#To call first:

#inter $bondId_drude drude $temp_core $gamma_core $temp_drude $gamma_drude $k_drude $mass_red_drude $r_cut

proc add_drude_to_core { bondId_drude id_core id_drude type_drude polarization } {

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

	set k_drude [lindex [inter $bondId_drude] 7]
	set mass_drude [lindex [inter $bondId_drude] 8]
	
	set q_core [part $id_core print q]
	set mass_core [part $id_core print mass]

	if {$q_core > 0} {
		set sign_q 1.0
	} else {
 		set sign_q -1.0
        }
#	set q_drude [expr $sign_q * pow($k_drude*$polarization/pow($sigma_core,3.0), 1./2.)]
	set q_drude [expr $sign_q * pow($k_drude*$polarization, 1./2.)]

	
	part $id_core q [expr $q_core - $q_drude] mass [expr $mass_core-$mass_drude] temp 0 gamma 0
	set drude_px [expr [lindex [part $id_core print pos] 0] + [t_random]*$polarization]
	set drude_py [expr [lindex [part $id_core print pos] 1] + [t_random]*$polarization]
	set drude_pz [expr [lindex [part $id_core print pos] 2] + [t_random]*$polarization]
	part $id_drude pos $drude_px $drude_py $drude_pz v 0 0 0  q $q_drude type $type_drude mass $mass_drude temp 0 gamma 0
	part $id_core bond $bondId_drude $id_drude
	#part $id_core bond $bondId_subt_elec_hdfkj4hrkjfy $id_drude

	return "Drude particle created with parameters:\nid: $id_drude\ncharge: $q_drude\nmass: $mass_drude\nCore parameters changed to:\nMass: [expr $mass_core-$mass_drude]\nCharge:[expr $q_core - $q_drude]\n$warnings"
}
