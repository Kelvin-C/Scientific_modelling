import raytracer as rt
import numpy as np

#Numbering the lenses
single_surface = 0
planoconvex = 1
bi_surface = 2

lens_type = 2 #type the lens's number to use it.

#Numbering the type of program to run.
single = 0
find_best_curvature = 1

run_type = 1 #type the program number to run it.

############# Ray Properties ##################

ray_centre_coordinate = [0, 0, 0]
ray_radius = 10
ray_number_of_circles = 7
ray_number_in_1st_circle = 10
ray_direction_vector    = [0, 0, 1]

###############################################

####### For a 1-time curvature test ###########

lens_curvature = -0.0
lens_2nd_curvature = -0.2 #for bisurface lens

###############################################

####### To find the best curvature  ###########

lens_curvature_min = 0.0001
lens_curvature_max = 0.1
lens_curvature_steps  = 20

#for bisurface lens
lens_2nd_curvature_min = -0.1
lens_2nd_curvature_max = -0.0001
lens_2nd_curvature_steps = 20

###############################################

######### General Lens's Properties ###########
#The aperture or depth may change to create the correct lens.
#Aperture takes priority over depth for bisurface lens.

lens_z0            = 10
lens_refractive_index   = 1.5168
surrounding_refractive_index = 1
lens_aperture   = 20    #Not required for Plano-convex Lens
lens_depth      = 20     #Not required for Spherical Lens

#The option to show change in aperture or depth. "Yes" or "No". By default, it is "No".
#If you are trying to find the best curvature, "Yes" will keep displaying the changes, and may slow the computer down.
show_change_in_depth_aperture = "Yes"   

#Methods of finding the focal length:
min_standard_deviation = 0      #Minimising standard deviation of focal point locations
max_number_intersections = 1    #Maximising the number of intersections at focal point.

#Type in the number that corresponds to the method for finding focal point.
focal_method = 0

##############  Functions  ####################

def propagateray(ray, lens_type, lens_curvature, lens_2nd_curvature = 0, word = "Yes"):
    """
    This function propagates the ray towards the lens. The RMS deviation at the focal point and the new ray object is returned.
    """
    if lens_type == 0:
        lens = rt.SphericalRefraction(lens_z0, lens_aperture, lens_curvature, surrounding_refractive_index, lens_refractive_index, word)    
    elif lens_type == 1:
        lens = rt.PlanoConvex(lens_z0, lens_depth, lens_curvature, surrounding_refractive_index, lens_refractive_index, word)
    elif lens_type == 2:
        lens = rt.BiSurface(lens_z0, lens_aperture, lens_depth, lens_curvature, lens_2nd_curvature, surrounding_refractive_index, lens_refractive_index, word)
    temp = ray.propagate(lens)
    if temp is None:
        return -1, -1
    if focal_method == 0:
        focal_position = ray.focalpoint_std()
    elif focal_method == 1:
        focal_position = ray.focalpoint_maxlength()
    if focal_position is None:
        output = rt.OutputPlane((lens_z0+lens_depth)*2.5)
        print "No focal point outside of the lens"
    else:
        output = rt.OutputPlane(focal_position[-1])
    temp = ray.propagate(output)
    if temp is None:
        return -1, -1
    rms = ray.rms_deviation("output")
    return rms, ray
    
def removerms(rmslist, raylist, curvaturelist):
    """
    Removes all of the RMS deviations that are 0.
    """
    rmslist, raylist, curvaturelist = removerms_helper(rmslist, raylist, curvaturelist, 0)
    return rmslist, raylist, curvaturelist

def removerms_helper(rmslist, raylist, curvaturelist, idx):
    """
    This function is intended to be ran from removerms
    """
    if len(rmslist) == 0:
        return rmslist, raylist, curvaturelist,    
    if idx < len(rmslist):
        if rmslist[idx] <= 0.000000001:
            del rmslist[idx]
            del raylist[idx]
            del curvaturelist[idx]
            removerms_helper(rmslist, raylist, curvaturelist, idx)
        else:        
            removerms_helper(rmslist, raylist, curvaturelist, idx+1)    
    return rmslist, raylist, curvaturelist
        
###############################################
    
############## Main Script ####################
    
run_type = int(run_type)
lens_type = int(lens_type)

if run_type == 0:
    ray_object = rt.RayArray(ray_centre_coordinate, ray_radius, ray_number_of_circles, ray_number_in_1st_circle, ray_direction_vector)
    if lens_type == 2:
        rms, ray = propagateray(ray_object, lens_type, lens_curvature, lens_2nd_curvature, show_change_in_depth_aperture)
    else:
        rms, ray = propagateray(ray_object, lens_type, lens_curvature, show_change_in_depth_aperture)
    if rms == -1:
        print "All of the rays were either total-interally reflected or they missed the lens"
    else:
        print "1st curvature =", lens_curvature
        if lens_type ==2:
            print "2nd curvature =", lens_2nd_curvature
        print "Ray Radius = ", ray_radius
        ray.plot()
        ray.spot()

elif run_type == 1:
    lens_curvature_array = np.linspace(lens_curvature_min, lens_curvature_max, lens_curvature_steps)
    lens_2nd_curvature_array = np.linspace(lens_2nd_curvature_min, lens_2nd_curvature_max, lens_2nd_curvature_steps)
    rmslist = []
    raylist = []
    curvaturelist = []
    if lens_type == 2:
        for curvature1 in lens_curvature_array:
            for curvature2 in lens_2nd_curvature_array:
                print "1st curvature =", curvature1
                print "2nd curvature =", curvature2
                ray_object = rt.RayArray(ray_centre_coordinate, ray_radius, ray_number_of_circles, ray_number_in_1st_circle, ray_direction_vector)
                rms, ray = propagateray(ray_object, lens_type, curvature1, curvature2, show_change_in_depth_aperture)
                    
                if rms != -1:
                    rmslist += rms,
                    raylist += ray,
                    curvaturelist += [[curvature1, curvature2]]
    else:
        for curvature in lens_curvature_array:
            print "curvature =", curvature
            ray_object = rt.RayArray(ray_centre_coordinate, ray_radius, ray_number_of_circles, ray_number_in_1st_circle, ray_direction_vector)
            rms, ray = propagateray(ray_object, lens_type, curvature, show_change_in_depth_aperture)
            if rms != -1:
                rmslist += rms,
                raylist += ray,
                curvaturelist += curvature,
                
    rmslist, raylist, curvaturelist = removerms(rmslist, raylist, curvaturelist)
    for i in range(len(rmslist)):
        if rmslist[i] == 0:
            rmslist.remove(rms)
    if rmslist == []:
        print "No curvature in the range were suitable"
    else:
        rmsindex = rmslist.index(min(rmslist))
        ray = raylist[rmsindex]        
        if lens_type ==2:
            print "Best 1st curvature =", curvaturelist[rmsindex][0]
            print "Best 2nd curvature =", curvaturelist[rmsindex][1]
        else:
            print "Best curvature =", curvaturelist[rmsindex]
        print "Ray Radius = ", ray_radius
        ray.plot()
        ray.spot()
