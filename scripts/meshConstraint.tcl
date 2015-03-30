
proc meshToParticles {path_to_mesh first_index ptype} {

	set infile [open "$path_to_mesh" "r"]
	set file_data [read $infile]
	close $infile

	set data [split $file_data "\n"]
	set mesh_coords ""
	set mesh_normals ""
	set tc 0
	set co ""
	foreach line $data {
		#TODO: COLLECT POSITIONS
		set inormal [string first "facet normal" $line] 
		if {$inormal != -1} {
			lappend mesh_normals [join [string range $line [expr $inormal + 12] end] " "] 
		} else {
			set ivertex [string first "vertex" $line] 
			if {$ivertex != -1} {
				lappend co [join [string range $line [expr $ivertex + 6] end] " "] 
				incr tc
				if {$tc==3} {
					set xc [expr ([lindex [lindex $co 0] 0] + [lindex [lindex $co 1] 0] + [lindex [lindex $co 2] 0])/3.0]
					set yc [expr ([lindex [lindex $co 0] 1] + [lindex [lindex $co 1] 1] + [lindex [lindex $co 2] 1])/3.0]
					set zc [expr ([lindex [lindex $co 0] 2] + [lindex [lindex $co 1] 2] + [lindex [lindex $co 2] 2])/3.0]
					lappend mesh_coords [list $xc $yc $zc] 
					set co ""
					set tc 0
				}
			}
		}
	}
	set meshParticles [llength $mesh_coords]

	for {set i 0} {$i < $meshParticles} {incr i} {
		set c [lindex $mesh_coords $i]
		part [expr $first_index+$i] pos [lindex $c 0] [lindex $c 1] [lindex $c 2] type $ptype fix 1 1 1
	}
	
	return "Created $meshParticles particles from file $path_to_mesh"
}



proc mesh_capacitor {list_path_to_mesh list_potentials list_types} { 
	
	set icc_areas [ list ]
	set icc_normals [ list ]
	set icc_epsilons [ list ]
	set icc_sigmas [ list ]
	set icclist ""

	set icc_area 1.0
	set icc_sigma 0
	set icc_eps 100000


	for {set i 0} {$i < $iccParticles} {incr i} {
		lappend icc_normals  [ list $nx $ny $nz ]
		lappend icc_areas $icc_area
		lappend icc_sigmas $icc_sigma 
		lappend icc_epsilons $icc_eps

		lappend icclist $i 
		if {[t_random] > 0.5} { set pm 1.0 } else { set pm -1.0}
		part $i type $icc_wall_type mass $m_carbon fix 1 1 1
	}

	return $res
}
