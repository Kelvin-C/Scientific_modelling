"""
Ray Tracer Module
"""
import numpy as np
import matplotlib.pyplot as plt

class Ray:
    """
    Creates an object which represents a ray of light.
    Parameters:
    position = 3-dimensional list
    direction = 3-dimensional list
    """
    
    def __init__(self, position = [0, 0, 0], direction = [1, 1, 1]):
        if len(position) != 3 or len(direction) != 3:
            raise Exception("This is a 3D Ray Tracer. It needs 3 dimensions")
        self.__position = np.array([position]).astype(float)
        self.__direction = np.array([direction]).astype(float)
        
    def __repr__(self):
        return """
        p = [%g, %g, %g]\n
        k = [%g, %g, %g]
        """ %(self.__position[-1][0], self.__position[-1][1], self.__position[-1][2], self.__direction[-1][0], self.__direction[-1][1], self.__direction[-1][2])
        
    def __str__(self):
        return """Current position of ray:\n
        position = [%g, %g, %g]\n
        direction = [%g, %g, %g]
        """ %(self.__position[-1][0], self.__position[-1][1], self.__position[-1][2], self.__direction[-1][0], self.__direction[-1][1], self.__direction[-1][2])
        
    def p(self, position = None):
        """
        Shows the current position of the ray.
        """
        return self.__position[-1]
        
    def k(self):
        """
        Shows the direction/velocity of the ray.
        """
        return self.__direction[-1]
        
    def append(self, position, direction):
        """
        Changes the direction/position of the ray.
        """
        self.__position = np.vstack([self.__position, position])
        self.__direction = np.vstack([self.__direction, direction])
               
    def _spot(self, rms_deviation):
        """
        Intended to be ran from the spot function of the RayArray class.
        This function shows the location of the ray on the x-y graph.
        For the parameter: type "input" to show the location at input, and "output" to show
        the location at output.
        """
        if rms_deviation == "input":
            plt.plot(self.__position[0][0], self.__position[0][1], 'bo')
        
        elif rms_deviation == "output":
            plt.plot(self.__position[-1][0], self.__position[-1][1], 'ro')
           
    def trace(self):
        """
        Traces the ray and plots it on a z-x diagram.
        """
        plt.plot(self.__position.T[2], self.__position.T[0],'b-')
        plt.show()
        
    def vertices(self):
        """
        Returns all locations of the rays.
        """
        return self.__position
        
class RayArray:
    """
    Generates an object which represents set of rays in a circular formation.
    The class creates a formation which contain circles expanding from the centre.
    Parameters:
    centre = the coordinates of the centre of the input rays (i.e. [0, 0, 0])
    radius = the maximum radius of the set of rays. (integer or float)
    no_of_circles = the number of circles expanding from the centre (integer)
    no_in_1st_circle = the number of rays which make up the first circle (integer)
    direction = the direction of each ray (assume all rays are parallel) (i.e. [0, 0, 0])
    """
    
    def __init__(self, centre = [0, 0, 0], radius = 10, no_of_circles = 7, no_in_1st_circle = 10, direction = [0, 0, 1]):
        self.__centre = centre
        self.__radius = float(radius)
        self.__no_of_circles =  no_of_circles
        self.__no_in_1st_circle = no_in_1st_circle
        self.__direction = direction
        self.__show_iter_number = 0
        
        distance_per_circle = float(radius)/no_of_circles
        ray_array=[Ray([centre[0], centre[1], centre[2]], direction)]
        for i in range(1,no_of_circles+1):
            no_rays_in_circle = no_in_1st_circle * i
            angle_per_ray = (2*np.pi)/(no_rays_in_circle)
            counter = 1
            while counter <= no_rays_in_circle:
                counter += 1
                ray_array += [Ray([centre[0]+i*distance_per_circle*np.cos(counter*angle_per_ray),centre[1]+i*distance_per_circle*np.sin(counter*angle_per_ray),centre[2]], direction)]
        self.__ray_array = ray_array
                                                                  
    def plot(self):
        """
        Traces all of the rays onto an z,x diagram.
        """
        plt.figure(figsize=(8,5))
        for i in self.__ray_array:
            Ray.trace(i)
        plt.title("Ray trace")
        plt.xlabel("z-axis")
        plt.ylabel("x-axis")
        plt.show()
              
    def propagate(self, surface):
        """
        This propagates the ray towards the chosen surface.
        If the ray does not hit the lens or it experiences total internal
        reflection, the ray is deleted.
        If this function returns None, then all of the rays are deleted. 
        """
        temp = self.__propagate_helper(surface, 0)
        return temp
               
    def __propagate_helper(self, surface, idx):
        if len(self.__ray_array) == 0:
            return None
        ray = self.__ray_array[idx]
        intersect = surface._propagate_ray(ray)
        if intersect is None:
            self.__ray_array.remove(ray)
            if idx < len(self.__ray_array):
                self.__propagate_helper(surface, idx)
        elif idx < len(self.__ray_array)-1:
            self.__propagate_helper(surface, idx+1)
        return 1

    def spot(self):
        """
        Creates a spot plot of the x,y locations of the rays.
        It generates 2 plots: one at input and one at output.
        """
        plt.figure(figsize=(10,10))
        rms_deviation_input = RayArray.rms_deviation(self, "input")
        rms_deviation_output = RayArray.rms_deviation(self, "output")
        
        plt.subplot(2,1,1)
        for ray in self.__ray_array:
            ray._spot("input")
        plt.title("z=0; RMS deviation = %g" %(rms_deviation_input))
        plt.xlabel("x-axis")
        plt.ylabel("y-axis") 
           
        plt.subplot(2,1,2) 
        for ray in self.__ray_array:
            ray._spot("output")
        z = self.__ray_array[0].p()[-1]
        plt.title("z=%g, focal point; RMS deviation = %g" %(z, rms_deviation_output))
        plt.xlabel("x-axis")
        plt.ylabel("y-axis") 
               
        plt.show()        
    
    def rms_deviation(self, place = "output"):
        """
        Finds the root-mean-squared deviation from the centre of the rays
        Parameter:
        place = type "input" for rms at input, or "output" for rms at output.
        """
        rms_list = []
        counter = 0.
        if place == "output":
            centrex = self.__ray_array[0].vertices()[-1][0]
            centrey = self.__ray_array[0].vertices()[-1][1]
            for i in range(1, len(self.__ray_array)):
                rms_list += (centrex - self.__ray_array[i].vertices()[-1][0])**2 + (centrey - self.__ray_array[i].vertices()[-1][1])**2,
                counter += 1
            rms = sum(rms_list)
        elif place == "input":
            centrex = self.__ray_array[0].vertices()[0][0]
            centrey = self.__ray_array[0].vertices()[0][1]
            for i in range(1, len(self.__ray_array)):
                rms_list += (centrex - self.__ray_array[i].vertices()[0][0])**2 + (centrey - self.__ray_array[i].vertices()[0][1])**2,
                counter += 1
            rms = sum(rms_list)
        if counter == 0:
            rms = 0
        else:
            rms = np.sqrt(rms/counter)
        return rms
                           
    def focalpoint_std(self):
        """
        Finds the focal point of the rays. It is solved by solving the vector equations
        and finding the directional vector constant. 
        Each ray's location and direction is compared to others to find the constant, and the focal
        point is found by minimising standard deviation of the constants.
        """
        std_list = []
        average_constant = []
        for i in range(len(self.__ray_array)):
            constant = []
            counter = 0.
            for j in range(len(self.__ray_array)):
                if j != i:
                    intersect_constant = self.intersect(i,j)
                    if intersect_constant != 12073 and intersect_constant > 0:
                        constant += intersect_constant,
                        counter += 1.
            if counter > 0:
                constant_array = np.array(constant)
                if np.std(constant) > 0.0000005:
                    std_list += [[np.std(constant), i]]
                average_constant += np.mean(constant_array),
        std_array = np.array(std_list)
        std_array_trans = std_array.T
        focal_position = [0,0,0]
        if len(std_list) == 0:
            return None
        else:
            std = min(std_array_trans[0])
            stdidx = int(np.where(std_array_trans[0] == std)[0][0])
            rayidx = int(std_array_trans[1][stdidx])
            focal_position = self.__ray_array[rayidx].p() + average_constant[stdidx]*self.__ray_array[rayidx].k()
            return focal_position
            
    def focalpoint_maxlength(self):
        """
        Finds the focal point of the rays. It is solved by solving the vector equations
        and finding the directional vector constant. 
        Each ray's location and direction is compared to others to find the constant, and the focal
        point is found by maximising the number of intersections between rays.
        """
        average_constant = []
        constant_length = []
        for i in range(len(self.__ray_array)):
            constant = []
            for j in range(len(self.__ray_array)):
                if j != i:
                    intersect_constant = self.intersect(i,j)
                    if intersect_constant != 12073 and intersect_constant > 0:
                        constant += intersect_constant,
            if len(constant) > 0:
                constant_length = [[len(constant), i]]
                average_constant += np.mean(constant),
        constant_length_array = np.array(constant_length)
        constant_length_trans = constant_length_array.T
        focal_position = [0,0,0]
        if len(constant_length) == 0:
            return None
        else:
            constant_length = max(constant_length_trans[0])
            constantidx = int(np.where(constant_length_trans[0] == constant_length)[0][0])
            rayidx = int(constant_length_trans[1][constantidx])
            focal_position = self.__ray_array[rayidx].p() + average_constant[constantidx]*self.__ray_array[rayidx].k()
            return focal_position
        
    def intersect(self, i, j):
        """
        Checks if two rays intersect and returns a constant, which should be multiplied by
        the 1st ray's direction vector to find the focal point position vector. 
        Parameters:
        i = The ray index of the 1st ray
        j = The ray index of the 2nd ray
        """
        
        p1 = self.__ray_array[j].p()  
        k1 = self.__ray_array[j].k()
        p2 = self.__ray_array[i].p()
        k2 = self.__ray_array[i].k()
        if k1[1] != 0:
            if (k2[0] - (k2[1]*k1[0]/k1[1])) != 0:
                    constant = ((k1[0]/k1[1])*(p2[1]-p1[1]) + p1[0] - p2[0])/(k2[0] - (k2[1]*k1[0]/k1[1]))
                    return constant
        return 12073 #Just a random number that's extremely unlikly to equal the constant.
        
class OpticalElement:
    
    def __init__(self, z0 = 100, n1 = 1, n2= 1.5):
        """
        Base object for lenses.
        """
        
        self.__z0 = float(z0)
        self.__n1 = float(n1)
        self.__n2 = float(n2)
    
    def n1(self):
        return self.__n1
            
    def n2(self):
        return self.__n2
        
    def z0(self):
        return self.__z0     
                          
    def _refraction(self, ray, normal, n1, n2):
        """
        This is intended to be ran by the _propagate_ray function from SphericalRefraction class.
        Refracts the ray when it iteracts with an optical element.
        This uses Snell's Law and linear algebra.
        Parameters:
        ray = The incident ray object.
        normal = The normal to the surface. The normal vector (3 dimensional list) should point into the surface from the incident ray side
        n1 = The refractive index of the medium in the incident ray side
        n2 = The refractive index of the medium in the refracted ray side
        """
        normal = normal.astype(float)
        if normal[-1] < 0:
            normal = -normal
        normal_hat = normal/np.linalg.norm(normal)
        k1_hat = ray.k()/np.linalg.norm(ray.k())
        k1mag = np.linalg.norm(ray.k())
        theta1 = np.arccos(np.dot(normal_hat,k1_hat))
        n_ratio = n1/n2
        
        if np.sin(theta1) > n2/n1:  #check for total-internal reflection
            return None
            
        theta2 = np.arcsin((n_ratio)*np.sin(theta1))
  
        if n1 < n2:
            check = (n_ratio**2+(-2*k1mag*n_ratio*np.cos(theta2)+1)/k1mag**2)
            if check < 0:
                return None
            else:
                k1_factor = np.sqrt(n_ratio**2+(-2*k1mag*n_ratio*np.cos(theta2)+1)/k1mag**2)
                k2 = k1_factor*ray.k() + normal_hat
            
        elif n1 > n2:
            check = (n_ratio**2+(-2*k1mag*n_ratio*np.cos(np.pi-theta2)+1)/k1mag**2)
            if check < 0:
                return None
            else:
                k1_factor = np.sqrt(n_ratio**2+(-2*k1mag*n_ratio*np.cos(np.pi-theta2)+1)/k1mag**2)
                k2 = k1_factor*ray.k() - normal_hat
            
        elif n1 == n2:
            return ray.k()
            
        return k2
        
        
class SphericalRefraction(OpticalElement):
    """
    Creates an optical element with a spherical or flat shape. This can be used to create
    a surface of the lens.
    Parameters:
    z0 = The z-coordinate of the centre of the surface (integer or float)
    aperture = The radial distance of the edge of the surface from the z-axis. (integer or float)
    curvature = The curvature of the surface (integer or float)
    n1 = The refractive index of the medium in the incident ray side (integer or float)
    n2 = The refractive index of the medium in the refracted ray side (integer or float)
    """
    
    def __init__(self, z0 = 20, aperture = 0.9, curvature = 1, n1 = 1, n2 = 1.5, word = "Yes"):
        OpticalElement.__init__(self, z0, n1, n2)
        self.__curvature = float(curvature)
        if curvature == 0:
            self.__radius = 0
            self.__centre = z0
            self.__aperture = float(aperture)
        else:
            self.__radius = 1./abs(curvature)
            if aperture >= self.__radius:
                self.__aperture = self.__radius - 0.0001
                if word == "Yes":
                    print "Aperture has been changed to %0.5f because chosen depth is too large" %self.__aperture
            else:
                self.__aperture = float(aperture)
            self.__centre = np.array([0, 0, z0+(1./self.__curvature)])
    
    def __repr__(self):
        return """z0 = %0.2f\n
        aperture =  %0.2f\n
        curvature =  %0.2f\n
        n1 =  %0.2f\n
        n2 =  %0.2f\n
        radius =  %0.2f\n""" %(self.z0(), self.aperture(), self.curvature(), self.n1(), self.n2(), self.radius())
        
    def __str__(self):
        return """These are the properties of the surface:\n
        The z-coordinate of the centre is at %0.2f;\n
        The aperture is %0.2f;\n
        The curvature is %0.2f;\n
        The refractive index of the medium on the incident ray side is %0.2f;\n
        The refractive index of the medium on the refracted ray side is %0.2f;\n
        The radius of curvature is %0.2f""" %(self.z0(), self.__aperture, self.__curvature, self.__n1(), self.n2(), self.__radius)
                   
    def __interceptcheck(self,intercept_position,ray):
        """
        Checks if the location of intercept is on the surface.
        Parameters:
        intercept_position = The 3-dimensional coordinate of the location of intercept.
        ray = The ray object
        """
        if self.__curvature == 0:
            distance_from_opticalaxis = np.sqrt(intercept_position[0]**2+intercept_position[1]**2)
            if distance_from_opticalaxis > self.__aperture:
                return None
            else:
                return intercept_position
                
        if list(intercept_position) == list(ray.p()):#if intercept position is at same place as current position
            return None                              #then it will ignore it.
            
        else: #This ensures the positions are actually on the surface.
            if self.__curvature < 0:
                radius_boundary = self.__centre[-1] + np.sqrt(self.__radius**2 - self.__aperture**2)
                if intercept_position[-1] < radius_boundary:
                    return None
                else:
                    return intercept_position
            else:
                radius_boundary = self.__centre[-1] - np.sqrt(self.__radius**2 - self.__aperture**2)
                if intercept_position[-1] > radius_boundary:
                    return None
                else:
                    return intercept_position
    
    def _intercept(self,ray):
        """
        Intended to be ran by _propagate_ray function.
        Find the location on the surface which the ray hits.
        Parameters:
        ray = The ray object.
        """
        
        k_hat = ray.k()/np.linalg.norm(ray.k())
        if self.__radius == 0:     #flat surface
            l = (self.z0() - ray.p()[-1])/k_hat[-1]
            intercept_position = ray.p()+l*k_hat
            intercept_position = self.__interceptcheck(intercept_position,ray)
            return intercept_position
                                     
        r = ray.p()-self.__centre
        r_dot_k_hat = np.dot(r,k_hat)
        determinant = r_dot_k_hat**2-(np.linalg.norm(r)**2-self.__radius**2)
        
        if determinant < 0:     #no intercept
            return None
            
        elif determinant <= 0.000005:   #only one intercept (took into account errors)
            l = -r_dot_k_hat
            intercept_position = ray.p()+l*k_hat
            intercept_position = self.__interceptcheck(intercept_position,ray)
            return intercept_position
            
        elif determinant > 0:  #Intercepts twice
            l = np.array([-r_dot_k_hat-(determinant)**0.5, -r_dot_k_hat+(determinant)**0.5])
            intercept_positions = np.array([ray.p()+l[0]*k_hat, ray.p()+l[1]*k_hat])
            true_intercept_positions = []
            for position in intercept_positions:
                temp = self.__interceptcheck(position,ray)
                if  temp is not None:
                    true_intercept_positions += [temp]
                    
            for intercept_position in true_intercept_positions:
                if intercept_position is not None:
                    return intercept_position
        
    def _propagate_ray(self, ray):
        """
        It checks if the ray hits the surface. If the ray hits, it will refract and change the
        direction of the ray. If it does not hit, the ray will be deleted.
        Parameters:
        ray = The ray object
        """
        intercept_position = self._intercept(ray)
        if intercept_position is not None:
            if self.__curvature == 0:
                normal = np.array([0,0,1])
            else:
                normal = self.__centre-intercept_position
            refracted_direction = self._refraction(ray, normal, self.n1(), self.n2())
            if refracted_direction is None:
                return None
            else:
                ray.append(intercept_position, refracted_direction)
                return 1
        return None

class PlanoConvex(SphericalRefraction):
    """
    Creates a plano-convex lens object. It uses the algorithms in SphericalRefraction
    to propagate rays.
    Parameters:
    z0 = The location of the surface, which the ray hits first, on the z-axis (integer or float)
    depth = The thickness of the lens (integer or float)
    curvature = The curvature of the lens (integer or float)
    n1 = The refractive index of the surroundings (integer or float)
    n2 = The refractive index of the lens (integer or float)
    """
    
    def __init__(self, z0 = 10, depth = 500, curvature = 1, n1 = 1, n2= 1.5, word = "Yes"):
        OpticalElement.__init__(self, z0, n1, n2) 
        self.__curvature = curvature
        if curvature == 0 or depth == 0:
            self.__aperture = 9999
            self.__flat_surface =  SphericalRefraction(z0, self.__aperture, curvature , n1, n2)
        else:
            radius = 1./abs(curvature)
            self.__radius = radius
            if depth >= radius:
                self.__depth = radius - 0.000001
                if word == "Yes":
                    print "Depth has been changed to %0.5f because chosen depth is too large" %self.__depth
            else:
                self.__depth = float(depth)
            self.__aperture = np.sqrt(radius**2-(self.__depth-radius)**2)
            if curvature < 0:
                self.__flat_surface =  SphericalRefraction(z0-self.__depth, self.__aperture, 0, n1, n2)
                self.__curved_surface = SphericalRefraction(z0, self.__aperture, curvature, n2, n1)
            else:
                self.__flat_surface =  SphericalRefraction(z0+self.__depth, self.__aperture, 0, n2, n1)
                self.__curved_surface = SphericalRefraction(z0, self.__aperture, curvature, n1, n2)
        SphericalRefraction.__init__(self, z0, self.__aperture, curvature, n1, n2)                      

    def __repr__(self):
        return """z0 = %0.2f\n
        depth =  %0.2f\n
        curvature =  %0.2f\n
        n1 =  %0.2f\n
        n2 =  %0.2f\n
        radius = %0.2f\n""" %(self.z0(), self.__depth, self.__curvature, self.n1(), self.n2(), self.__radius)
            
    def __str__(self):
        return """These are the properties of the plano-convex lens:\n
        The z-coordinate of the surface closest to [0,0,0] is at %0.2f;\n
        The aperture is %0.2f;\n
        The curvature is %0.2f;\n
        The refractive index of the surroundings is %0.2f;\n
        The refractive index of the lens %0.2f;\n
        The radius of curvature is %0.2f""" %(self.z0(), self.__aperture, self.__curvature, self.__n1(), self.n2(), self.__radius)
            
    def _propagate_ray(self, ray):
        """
        It checks if the ray hits the surface. If the ray hits, it will refract and change the
        direction of the ray. If it does not hit, the ray will be deleted.
        Parameters:
        ray = The ray object
        """
        if self.__curvature == 0:
            a = SphericalRefraction._propagate_ray(self.__flat_surface, ray)
            b = 5
        elif self.__curvature > 0:
            a = SphericalRefraction._propagate_ray(self.__curved_surface, ray)
            b = SphericalRefraction._propagate_ray(self.__flat_surface, ray)
        elif self.__curvature < 0:
            a = SphericalRefraction._propagate_ray(self.__flat_surface, ray)           
            b = SphericalRefraction._propagate_ray(self.__curved_surface, ray)
        if a is None or b is None:
            return None
        else:
            return 1
            
class BiSurface(PlanoConvex):
    """
    This creates an object that simulates a lens with a flat or curved surface on both sides.
    To make the object, you need:
    z0 = the z-coordinate of the centre of the surface that is closest to [0,0,0]. 
    aperture = the radial distance of the surface from the z-axis (optical axis).
    depth = the thickness of the lens.
    curvature1 = the curvature of the surface closest to [0,0,0] (+ve curvature means the surface is curved towards -z.)
    curvature2 = the curvature of the opposite surface. (+ve curvature means the surface is curved towards -z.)
    n1 = surrounding refractive index
    n2 = lens refractive index
    
    The aperture and depth will automatically be adjusted if they are unsuitable.
    When needed (i.e. biconvex lens), aperture takes priority over depth.
    """
    
    def __init__(self, z0 = 10, aperture = 500, depth = 2, curvature1 = 1, curvature2 = 1, n1 = 1, n2= 1.5, word = "Yes"):
            OpticalElement.__init__(self, z0, n1, n2)       
            
            if curvature1 == 0 and curvature2 == 0:
                self.__radius1 = "Flat Surface"
                self.__radius2 = "Flat Surface"
                                    
            elif curvature1 == 0:     #Situations for a flat surface as the 1st lens
                #No depth requirements needed for curvature2 < 0
                radius2 = 1./abs(curvature2)
                self.__radius1 = "Flat Surface"
                self.__radius2 = radius2
                if aperture >= radius2: #prevents the surface from making a full circle
                    aperture = radius2 - 0.00001
                    if word == "Yes":
                        print "Aperture has been changed to %0.5f because chosen aperture is too large" %aperture
                if curvature2 < 0:
                    depth = radius2 - np.sqrt(radius2**2 - aperture**2)
                    if word == "Yes":
                        print "Depth has been changed to %0.5f to match the aperture" %depth
                
            elif curvature2 == 0:   #Situations for a flat surface as the 2nd lens
                #No depth requirements needed for curvature1 < 0
                radius1 = 1./abs(curvature1)
                self.__radius2 = "Flat Surface"
                self.__radius1 = radius1
                if aperture >= radius1: #prevents the surface from making a full circle
                    aperture = radius1 - 0.00001
                    if word == "Yes":
                        print "Aperture has been changed to %0.5f because chosen aperture is too large" %aperture
                if curvature1 < 0:
                    depth = radius1 - np.sqrt(radius1**2 - aperture**2)  
                    if word == "Yes":
                        print "Depth has been changed to %0.5g to match the aperture" %depth        
                
            elif curvature1 > 0 or curvature1 < 0:
                #Only when curvature2 < 0 is when the surface has a maxiumum depth.
                radius1 = 1./abs(curvature1)
                radius2 = 1./abs(curvature2)
                self.__radius1 = radius1
                self.__radius2 = radius2
                if aperture >= radius1: #prevents the lens from making a full circle
                    aperture = radius1 - 0.00001
                    if word == "Yes":
                        print "Aperture has been changed to %0.5f because chosen aperture is too large" %aperture
                if aperture >= radius2: #prevents the lens from making a full circle
                    aperture = radius2 - 0.00001
                    if word == "Yes":
                        print "Aperture has been changed to %0.5f because chosen aperture is too large" %aperture
                if curvature1 > 0 and curvature2 < 0:
                    depth = radius2 + radius1 - np.sqrt(radius2**2 - aperture**2) - np.sqrt(radius1**2 - aperture**2)
                    if word == "Yes":
                        print "Depth has been changed to %0.5g to match the aperture" %depth
                
            self.__surface1 = SphericalRefraction(z0, aperture, curvature1 , n1, n2)
            self.__surface2 = SphericalRefraction(z0+depth, aperture, curvature1 , n2, n1) 
            self.__curvature1 = curvature1
            self.__curvature2 = curvature2     
            self.__depth = depth
            self.__aperture = aperture
   
    def __repr__(self):
        return """z0 = %0.2f\n
        aperture = %0.2f\n
        depth =  %0.2f\n
        curvature1 =  %0.2f\n
        curvature2 =  %0.2f\n
        n1 =  %0.2f\n
        n2 =  %0.2f\n
        radius1 = %0.2f\n
        radius2 = %0.2f""" %(self.z0(), self.__aperture, self.__depth, self.__curvature1, self.__curvature2, self.n1(), self.n2(), self.__radius1, self.__radius2)
            
    def __str__(self):
        return """These are the properties of the plano-convex lens:\n
        The z-coordinate of the surface closest to [0,0,0] is at %0.2f;\n
        The depth is %0.2f;\n
        The aperture is %0.2f;\n
        The curvature of the surface closest to [0,0,0] is %0.2f;\n
        The curvature of the surface furthest from [0,0,0] is %0.2f;\n
        The refractive index of the surroundings is %0.2f;\n
        The refractive index of the lens %0.2f;\n
        The radius of curvature of the surface closest to [0,0,0] is %0.2f;\n
        The radius of curvature of the furthest closest to [0,0,0] is %0.2f""" %(self.z0(), self.__aperture, self.__depth, self.__curvature1, self.__curvature2, self.n1(), self.n2(), self.__radius1, self.__radius2)
              
                  
                    
    def _propagate_ray(self, ray):
        """
        It checks if the ray hits the surface. If the ray hits, it will refract and change the
        direction of the ray. If it does not hit, the ray will be deleted.
        Parameters:
        ray = The ray object
        """
        a = SphericalRefraction._propagate_ray(self.__surface1, ray)
        b = SphericalRefraction._propagate_ray(self.__surface2, ray)
        if a is None or b is None:
            return None
        else:
            return 1
        
class OutputPlane:
    """
    Creates a plane parallel to the x-y plane which stops the rays. This acts as an output.
    Parameters:
    z = The location of the plane on the z-axis
    """
    
    def __init__(self, z = 200):
        self.__z = float(z)
    
    def __repr__(self):
        return "z = %g" %self.__z
        
    def __str__(self):
        return "The output object has z-coordinate of the plane is %g" %self.__z
        
    def _propagate_ray(self, ray):
        """
        It checks if the ray hits the surface. If the ray hits, it will refract and change the
        direction of the ray. If it does not hit, the ray will be deleted.
        Parameters:
        ray = The ray object
        """
        k_factor = (self.__z-ray.p()[-1])/ray.k()[-1]
        output_x = ray.p()[0]+k_factor*ray.k()[0]
        output_y = ray.p()[1]+k_factor*ray.k()[1]
        output_z = self.__z
        output_p = np.array([output_x, output_y, output_z])
        ray.append(output_p, ray.k())
        return 1
        
