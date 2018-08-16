———————————————————————————————————————

raytracer

———————————————————————————————————————

Open the script called "script.py" with notepad or another text editor.
When the script is ran, it will open the corresponding ray trace and the spot diagram.

The variables and values that can be changed are the ones above the 'functions'.
However, variables that are used to number the lenses, type of program, focal
point methods or aperture/depth change display should not be changed.

———————————————————————————————————————

Definitions:

#Numbering the lenses
'single_surface' lens have infinitesimal thickness. [DO NOT CHANGE]
'planoconvex' lens have one surface flat and the other is curved. [DO NOT CHANGE]
'bi_surface' lens is the most general lens. Its surfaces can be curved or flat. [DO NOT CHANGE]
'lens_type' chooses which lens you use. [CHOOSE 0, 1 OR 2]

#Numbering the type of program to run.
'single' only run traces the rays for a single curvature. [DO NOT CHANGE]
'find_best_curvature' is the program which attempts to find the best curvature for the chosen lens. [DO NOT CHANGE]
'run_type' chooses which program to run. [CHOOSE 0 OR 1]

############# Ray Properties ##################
The rays are formed in a circular formation. The formation consists
of circles which surround the central ray.

'ray_centre_coordinate' is the coordinate of the central ray. [3D LIST ONLY]
'ray_radius' is the radius of the set of the rays [FLOAT OR INTEGER]
'ray_number_of_circles' is a value which defines the number of circles in the formation [INTEGER ONLY]
'ray_number_in_1st_circle' is a value which defines the number of rays that makes up the first circle around the centre. [INTEGER ONLY]
'ray_direction_vector' is a direction vector which the rays start with. [3D LIST ONLY]

###############################################

####### For a 1-time curvature test ###########
These curvature settings are only used when the 'single' program is ran.

'lens_curvature' is the curvature of the surface closest to [0, 0, 0] coordinate. [INTEGER OR FLOAT]
'lens_2nd_curvature' is the curvature of the surface furthest from [0, 0, 0] coordinate. [INTEGER OR FLOAT]

###############################################

####### To find the best curvature  ###########
These curvature settings are only used when 'find_best_curvature' is ran.

'lens_curvature_min' is the minimum curvature of the lens closest to [0, 0, 0] coordinate. [INTEGER OR FLOAT]
'lens_curvature_max' is the maximum curvature of the lens closest to [0, 0, 0] coordinate. [INTEGER OR FLOAT]
'lens_curvature_steps' is a value which defines the number of curvatures between the range. [INTEGER ONLY]

#for bisurface lens 
'lens_2nd_curvature_min' is the minimum curvature of the lens furthest from [0, 0, 0] coordinate. [INTEGER OR FLOAT]
'lens_2nd_curvature_max' is the maximum curvature of the lens furthest from [0, 0, 0] coordinate. [INTEGER OR FLOAT]
'lens_2nd_curvature_steps' is a value which defines the number of curvatures between the range. [INTEGER ONLY]

###############################################

######### General Lens's Properties ###########
#The aperture or depth may change to create the correct lens.
#Aperture takes priority over depth for bisurface lens.

'lens_z0' is the z-coordinate of the centre of the lens closest to [0, 0, 0] coordinate [INTEGER OR FLOAT]
'lens_refractive_index' is the refractive index of the lens. [INTEGER OR FLOAT]
'surrounding_refractive_index' is the refractive index of the medium surrounding the lens. [INTEGER OR FLOAT]
'lens_aperture' is the aperture radius of the lens. [INTEGER OR FLOAT]
'lens_depth' is the thickness of the lens. [INTEGER OR FLOAT]

#The option to show change in aperture or depth. "Yes" or "No". By default, it is "Yes".
#If you are trying to find the best curvature, "Yes" will keep displaying the changes, and may slow the computer down.
'show_change_in_depth_aperture' is a setting which controls the display of the aperture and depth changes. ["Yes" or "No" ONLY] 

#Methods of finding the focal length:
'min_standard_deviation' is a method for finding the focal length by minimising the standard deviations of focal point positions.
'max_number_intersections' is a method for finding the focal length by maximising the number of intersections between rays at the point.
'focal_method' defines the method for finding the focal length. [CHOOSE 0 OR 1 ONLY]

———————————————————————————————————————

Contact details:

- kelvin.chan14@imperial.ac.uk
