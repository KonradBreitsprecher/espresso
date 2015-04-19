
proc meshToParticles {path_to_mesh first_index ptype} {
	global mesh_normals mesh_areas

	set infile [open "$path_to_mesh" "r"]
	set file_data [read $infile]
	close $infile

	set data [split $file_data "\n"]
	set mesh_coords ""
	set mesh_normals ""
	set mesh_areas ""
	set tc 0
	set co ""
	foreach line $data {
		set inormal [string first "facet normal" $line] 
		if {$inormal != -1} {
			lappend mesh_normals [join [string range $line [expr $inormal + 12] end] " "] 
		} else {
			set ivertex [string first "vertex" $line] 
			if {$ivertex != -1} {
				lappend co [join [string range $line [expr $ivertex + 6] end] " "] 
				incr tc
				if {$tc==3} {
					set A [lindex $co 0]
					set B [lindex $co 1]
					set C [lindex $co 2]

					set xc [expr ([lindex $A 0] + [lindex $B 0] + [lindex $C 0])/3.0]
					set yc [expr ([lindex $A 1] + [lindex $B 1] + [lindex $C 1])/3.0]
					set zc [expr ([lindex $A 2] + [lindex $B 2] + [lindex $C 2])/3.0]
					lappend mesh_coords [list $xc $yc $zc]
					
					lappend mesh_areas [expr 0.5 * [veclen [veccross_product3d [vecsub $B $A] [vecsub $C $A]]]] 
				
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
	
	return $meshParticles
}



proc mesh_capacitor_icc {firstIndex list_pathToMesh list_potentials list_types bins ext_pot_path} { 
	global mesh_normals mesh_areas
	global icc_areas icc_normals icc_epsilons icc_sigmas

	#CREATE PARTICLE FROM STL-MESHFILES AND SETUP ICC LISTS
	set icc_areas [ list ]
	set icc_normals [ list ]
	set icc_epsilons [ list ]
	set icc_sigmas [ list ]

	set icc_sigma 0
	set icc_eps 100000

	set num_particles [ list ]
	set i $firstIndex
	foreach pathToMesh $list_pathToMesh type $list_types {
		lappend num_particles [meshToParticles $pathToMesh $i $type] 
		foreach m $mesh_areas n $mesh_normals {
			lappend icc_areas $m
			lappend icc_normals $n
			lappend icc_sigmas $icc_sigma
			lappend icc_epsilons $icc_eps
			
			if {[t_random] > 0.5} { set pm 1.0 } else { set pm -1.0}
			part $i q [expr $pm*(0.001 + 0.001*[t_random])]
			incr i
		}
	}
	#CALCULATE EXTERNAL POTENTIAL AND SAVE TO FILE
	if {[file exists $ext_pot_path]} {
		puts "Found $ext_pot_path, skipping potential calculations"
	} else {
		puts [generate_potential_from_mesh [llength $list_pathToMesh] $list_pathToMesh $list_potentials $bins $ext_pot_path]
	    puts "Created $ext_pot_path" 
	}
	puts "Setup ICC lists, assigned small random charge +/- \[0.001:0.002\] to particles $firstIndex to [expr $i-1]\n"
	return $num_particles 
}
