"""Classes for modelling propellers."""

import os
import json

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import airfoil_db as adb

from abc import abstractmethod

from .poly_fit import poly_func
from .helpers import to_rpm, to_rads, import_value, check_filepath

class BaseProp:
    """Defines a propeller."""

    def __init__(self):
        self.name = None
        self.manufacturer = None
        self.diameter = None
        self.pitch = None

        
    def __str__(self):
        string = "Propeller: {0}".format(self.name)
        string += "\n\tManufacturer: {0}".format(self.manufacturer)
        string += "\n\tDiameter: {0} inches".format(self.diameter)
        string += "\n\tPitch: {0} inches".format(self.pitch)
        return string


    def get_advance_ratio(self, w, V):
        """Determines the advance ratio from angular velocity in rad/s and velocity in ft/s.

        Parameters
        ----------
        w : float
            Angular velocity in rad/s.

        V : float
            Velocity in ft/s.

        Returns
        -------
        float
            Advance ratio.
        """
        rps = to_rpm(w)/60.0
        if abs(rps) < 1e-10:
            return np.inf
        else:
            return 12.0*V/(rps*self.diameter)


    def get_velocity(self, rpm, J):
        """Determines the velocity from the rpm and advance ratio.

        Parameters
        ----------
        rpm : float
            Angular velocity in rpm.

        J : float
            Advance ratio.

        Returns
        -------
        float
            Velocity in ft/s.
        """
        rps = rpm/60.0
        return J*rps*self.diameter/12.0


    def plot_coefs(self, rpm_lims=[0.0, 35000.0], J_lims=[0.0, 1.4]):
        """Plot thrust and torque coefficients as functions of rpm and advance ratio.

        Parameters
        ----------
        rpm_lims : list, optional
            Limits for plotting in rpm. Defaults to [0.0, 35000.0].

        J_lims : list, optional
            Limits for plotting in advance ratio. Defaults to [0.0, 1.4].
        """

        # Get distributions
        rpms = np.linspace(rpm_lims[0], rpm_lims[1], 10)
        Js = np.linspace(J_lims[0], J_lims[1], 10)

        # Initialize figure
        fig = plt.figure(figsize=plt.figaspect(1.))
        fig.suptitle(self.name)
        ax = fig.add_subplot(1,2,1, projection='3d')

        # Loop through rpm values to plot thrust
        for rpm in rpms:
            thrust = np.zeros_like(Js)
            for i, J in enumerate(Js):
                w = to_rads(rpm)
                V = self.get_velocity(rpm, J)
                thrust[i] = self.get_thrust_coef(w, V)
            rpms_to_plot = np.full(len(thrust), rpm)
            ax.plot(Js, rpms_to_plot, thrust, 'r-')

        # Format figure
        ax.set_title("Predicted Thrust")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Thrust Coefficient")

        ax = fig.add_subplot(1,2,2, projection='3d')

        # Loop through rpm values to plot torque
        for rpm in rpms:
            power = np.zeros_like(Js)
            for i, J in enumerate(Js):
                w = to_rads(rpm)
                V = self.get_velocity(rpm, J)
                power[i] = self.get_torque_coef(w, V)*2.0*np.pi
            rpms_to_plot = np.full(len(power), rpm)
            ax.plot(Js, rpms_to_plot, power, 'r-')

        # Format
        ax.set_title("Predicted Power")
        ax.set_xlabel("Advance Ratio")
        ax.set_ylabel("RPM")
        ax.set_zlabel("Power Coefficient")
        plt.show()


    @abstractmethod
    def get_torque_coef(self, w, V):
        pass


    @abstractmethod
    def get_thrust_coef(self, w, V):
        pass


class DatabaseFitProp(BaseProp):
    """Defines a propeller by database fits.

    Parameters
    ----------
    record : list
        SQL database record defining this prop. This class should
        be initialized using the pyprop.create_component_from_database()
        function.

    """
    
    def __init__(self, record):

        # Store params
        self.name = record[1]
        self.manufacturer = record[2]
        self.diameter = float(record[3])
        self.pitch = float(record[4])

        # Determine fit properties
        thrust_fit_order = int(record[5])
        fit_of_thrust_fit_order = int(record[6])
        power_fit_order = int(record[7])
        fit_of_power_fit_order = int(record[8])
        num_coefs = (thrust_fit_order+1)*(fit_of_thrust_fit_order+1)
        self.thrust_coefs = record[9:num_coefs+9].reshape((thrust_fit_order+1,fit_of_thrust_fit_order+1)).astype(np.float)
        self.power_coefs = record[num_coefs+9:].reshape((power_fit_order+1,fit_of_power_fit_order+1)).astype(np.float)


    def get_torque_coef(self, w, V):
        """Returns the torque coefficient for the prop at the given angular velocity and freestream velocity.

        Parameters
        ----------
        w : float
            Angular velocity of the prop in radians per second.

        V : float
            Freestream velocity in feet per second.

        Returns
        -------
        float
            Torque coefficient.
        """

        # Get revs per minute and J
        rpm = to_rpm(w)
        J = self.get_advance_ratio(w, V)

        # Evaluate function to get fit coefficients
        a = poly_func(self.power_coefs.T, rpm)
        if(a[-1]>0): # Quadratic coefficient should always be non-positive
            a[-1] = 0

        # Evaluate function to get torque coefficient
        return poly_func(a, J)/2.0*np.pi
        

    def get_thrust_coef(self, w, V):
        """Returns the thrust coefficient for the prop at the given angular velocity and freestream velocity.

        Parameters
        ----------
        w : float
            Angular velocity of the prop in radians per second.

        V : float
            Freestream velocity in feet per second.

        Returns
        -------
        float
            Thrust coefficient.
        """

        # Get revs per minute and J
        rpm = to_rpm(w)
        J = self.get_advance_ratio(w, V)
        
        # Evaluate function to get fit coefficients
        a = poly_func(self.thrust_coefs.T, rpm)

        # Quadratic coefficient should always be non-positive
        if(a[-1]>0):
            a[-1] = 0

        return poly_func(a, J)


class DatabaseDataProp(BaseProp):
    """Defines a propeller by tabulated data. The data for these
    props are contained in the pyprop/props/ directory. PyProp will
    import that data automatically to intialize the prop.

    Parameters
    ----------
    name : str
        The name of the prop as stored in the props/ directory. This
        is the same as the name of the prop in the components database.
        The name should be given without the ".ppdat" or ".ppinf" file
        extension.
    """

    def __init__(self, name):
        
        # Load data file
        prop_dir = os.path.join(os.path.dirname(__file__), "props")
        data_filename = os.path.join(prop_dir, name+".ppdat")
        with open(data_filename, 'r') as data_file:
            self._data = np.genfromtxt(data_file, skip_header=1)

        # Load info
        info_filename = os.path.join(prop_dir, name+".ppinf")
        with open(info_filename, 'r') as info_file:
            self._info_dict = json.load(info_file)

        # Store params
        self.name = self._info_dict["name"]
        self.manufacturer = self._info_dict["manufacturer"]
        self.diameter = self._info_dict["diameter"]
        self.pitch = self._info_dict["pitch"]


    def get_torque_coef(self, w, V):
        """Returns the torque coefficient for the prop at the given angular velocity and freestream velocity.

        Parameters
        ----------
        w : float
            Angular velocity of the prop in radians per second.

        V : float
            Freestream velocity in feet per second.

        Returns
        -------
        float
            Torque coefficient.
        """

        # Get revs per minute and J
        rpm = to_rpm(w)
        J = self.get_advance_ratio(w, V)

        return interp.griddata(self._data[:,:2], self._data[:,3], np.array([[rpm, J]]), method='linear').item()/(2.0*np.pi)
        

    def get_thrust_coef(self, w, V):
        """Returns the thrust coefficient for the prop at the given angular velocity and freestream velocity.

        Parameters
        ----------
        w : float
            Angular velocity of the prop in radians per second.

        V : float
            Freestream velocity in feet per second.

        Returns
        -------
        float
            Thrust coefficient.
        """

        # Get revs per minute and J
        rpm = to_rpm(w)
        J = self.get_advance_ratio(w, V)
        
        return interp.griddata(self._data[:,:2], self._data[:,2], np.array([[rpm, J]]), method='linear').item()


class BladeElementProp(BaseProp):
    """Defines the performance of a propeller using blade element theory.

    Parameters
    ----------
    name : str
        Name of the prop.

    input_dict : dict
        Dictionary describing the geometry of the prop.

    airfoil_dict : dict, optional
        Dictionary of airfoils to be used by the prop. This is only needed
        if these have not already been given in input_dict.
    """

    def __init__(self, name, input_dict, airfoil_dict=None):

        # Initialize storage
        self._getter_data = {}

        # Load basic parameters
        self.name = name
        self._input_dict = input_dict
        geom_dict = self._input_dict["geometry"]
        self.diameter = import_value("diameter", geom_dict, "English", None)
        self.k = self._input_dict["geometry"]["n_blades"]
        self._r_hub = import_value("hub_radius", geom_dict, "English", 0.0)
        self.weight = import_value("weight", geom_dict, "English", 0.0)
        self._rot_dir = geom_dict.get("rotation", "CW")

        # Create distributions
        self._initialize_zeta()
        self._initialize_airfoils(airfoil_dict)
        self._initialize_getters()
        self._initialize_distributions()

    
    def _initialize_zeta(self):

        # Get grid parameters
        self._N = self._input_dict["geometry"]["grid"]
        distribution = self._input_dict["geometry"].get("distribution", "linear")

        # Determine nondimensional hub radius
        self._zeta_hub = 2.0*self._r_hub/self.diameter

        # Create array of nondimensional locations
        if distribution == "cosine_cluster": # Cosine clustering
            theta = np.linspace(0.0, np.pi, self._N)
            self._zeta = self._zeta_hub+0.5*(1-np.cos(theta))*(0.9999-self._zeta_hub)

        elif distribution == "linear": # Linear spacing
            self._zeta = np.linspace(self._zeta_hub, 0.9999, self._N)

        elif isinstance(distribution, list): # User-specified distribution
            if len(distribution) != self._N:
                raise IOError("User specified distribution must have length of N. Got length {0}; needed length {1}.".format(len(distribution), self._N))
            if abs(distribution[0]-self._zeta_hub)>1e-10:
                raise IOError("Distribution of zeta must begin at the hub. The value of zeta at the jub is {0}.".format(self._zeta_hub))

            self._zeta = np.array(distribution)

        else:
            raise IOError("Distribution type {0} not recognized for propeller {1}.".format(distribution, self.name))


    def _initialize_getters(self):
        # Sets getters for functions which are a function of span

        # Determine how the wing LQC has been given
        chord_data = import_value("chord", self._input_dict, "English", "not_given")
        pitch_data = import_value("pitch", self._input_dict, "English", "not_given")
        twist_data = import_value("twist", self._input_dict, "English", "not_given")
        selig_data = import_value("geometry_file", self._input_dict, "English", "not_given")

        # Check for redundant definitions
        if chord_data != "not_given" and isinstance(selig_data, np.ndarray):
            raise IOError("'chord' and 'geometry_file' may not both be specified at once.")
        if pitch_data != "not_given" and isinstance(selig_data, np.ndarray):
            raise IOError("'pitch' and 'geometry_file' may not both be specified at once.")
        if twist_data != "not_given" and isinstance(selig_data, np.ndarray):
            raise IOError("'twist' and 'geometry_file' may not both be specified at once.")
        if twist_data != "not_given" and pitch_data != "not_given":
            raise IOError("'twist' and 'pitch' may not both be specified at once.")

        # Selig data file
        if selig_data != "not_given":

            # Extract distributions
            raw_zeta = selig_data[:,0]
            raw_chord = selig_data[:,1]
            raw_twist = selig_data[:,2]

            # Package as if given by the user
            chord_data = np.concatenate(1, (raw_zeta, raw_chord))
            twist_data = np.concatenate(1, (raw_zeta, raw_twist))

        # Chord
        if chord_data == "not_given":
            raise IOError("'chord' is not optional if 'geometry_file' is not given.")
        if isinstance(chord_data, tuple): # Elliptic distribution
            self.get_chord = self._build_elliptic_chord_dist(chord_data[1])
        if callable(chord_data):
            self.get_chord = chord_data
        else:
            self.get_chord = self._build_getter_linear_f_of_span(chord_data, "chord", angular_data=False)

        # Set default twist of zero
        if pitch_data == "not_given" and twist_data == "not_given":
            twist_data = 0.0

        # Twist
        if twist_data != "not_given":
            if callable(twist_data):
                self.get_twist = twist_data
            else:
                self.get_twist = self._build_getter_linear_f_of_span(twist_data, "twist", angular_data=True)

        # Pitch
        if pitch_data != "not_given":
            if callable(pitch_data):
                self.get_pitch = pitch_data
            else:
                self.get_pitch = self._build_getter_linear_f_of_span(pitch_data, "pitch", angular_data=False)


    def _build_getter_linear_f_of_span(self, data, name, angular_data=False, flip_sign=False):
        # Defines a getter function for data which is a function of span

        if isinstance(data, float): # Constant
            if angular_data:
                self._getter_data[name] = m.radians(data)
            else:
                self._getter_data[name] = data

            def getter(span):
                """
                span : float or ndarray
                    Non-dimensional span location.
                """

                # Make input an array
                converted = False
                if isinstance(span, float):
                    converted = True
                    span = np.asarray(span)[np.newaxis]

                # Reverse sign
                if flip_sign:
                    data = -np.full(span.shape, self._getter_data[name])
                else:
                    data = np.full(span.shape, self._getter_data[name])

                # Convert back to scalar if needed
                if converted:
                    span = span.item()
                    return data.item()
                else:
                    return data

        
        else: # Array
            if isinstance(data[0], np.void): # This will happen if the user inputs the array params as ints
                new_data = np.zeros((data.shape[0],2), dtype=float)
                for i in range(data.shape[0]):
                    new_data[i,0] = data[i][0]
                    new_data[i,1] = data[i][1]
                data = new_data

            self._getter_data[name] = np.copy(data)

            def getter(span):
                """
                span : float or ndarray
                    Non-dimensional span location.
                """

                # Convert input to array
                converted = False
                if isinstance(span, float):
                    converted = True
                    span = np.asarray(span)[np.newaxis]

                # Perform interpolation
                if angular_data:
                    data = np.interp(span, self._getter_data[name][:,0], np.radians(self._getter_data[name][:,1]))
                else:
                    data = np.interp(span, self._getter_data[name][:,0], self._getter_data[name][:,1])

                # Reverse data
                if flip_sign:
                    data = -data

                # Convert back to scalar if needed
                if converted:
                    span = span.item()
                    return data.item()
                else:
                    return data

        return getter


    def _build_elliptic_chord_dist(self, root_chord):
        # Creates a getter which will return the chord length as a function of span fraction according to an elliptic distribution
        self._root_chord = root_chord

        def getter(span_frac):
            return self._root_chord*np.sqrt(1-span_frac*span_frac)

        return getter


    def _initialize_airfoils(self, airfoil_input):
        # Creates a dictionary of all the airfoils. This dictionary is then used to set up
        # the airfoil distribution.

        # Get airfoils from somewhere
        airfoil_input = self._input_dict.get("airfoils", airfoil_input)
        if airfoil_input is None:
            airfoil_input = { "default" : {} }

        # Load airfoil database from separate file
        if isinstance(airfoil_input, str):
            check_filepath(airfoil_input, ".json")
            with open(airfoil_input, 'r') as airfoil_db_handle:
                airfoil_input = json.load(airfoil_db_handle)

        # Make sure they gave us a dictionary
        if not isinstance(airfoil_input, dict):
            raise IOError("'airfoils' must be a string or dict.")

        # Get which airfoils are specified for this prop
        default_airfoil = list(airfoil_input.keys())[0]
        airfoil = import_value("airfoil", self._input_dict, "English", default_airfoil)

        self._airfoils = []
        self._airfoil_spans = []
        self._num_airfoils = 0

        # Setup data table
        if isinstance(airfoil, str): # Constant airfoil

            if not airfoil in list(airfoil_input.keys()):
                raise IOError("'{0}' must be specified in 'airfoils'.".format(airfoil))

            # Just put the same airfoil at the root and the tip
            self._airfoils.append(adb.Airfoil(airfoil, airfoil_input[airfoil]))
            self._num_airfoils = 1


        elif isinstance(airfoil, np.ndarray): # Distribution of airfoils

            # Store each airfoil and its span location
            for row in airfoil:

                name = row[1].item()

                try:
                    self._airfoils.append(adb.Airfoil(airfoil, airfoil_input[airfoil]))
                except NameError:
                    raise IOError("'{0}' must be specified in 'airfoils'.".format(name))

                self._airfoil_spans.append(float(row[0]))
                self._num_airfoils += 1

        else:
            raise IOError("Airfoil definition must a be a string or an array.")

        self._airfoil_spans = np.asarray(self._airfoil_spans)


    def _initialize_distributions(self):
        # Creates and stores vectors of important data at each control point

        self.u_a_cp = self._get_axial_vec(self._zeta)
        self.u_n_cp = self._get_normal_vec(self._zeta)
        self.u_s_cp = self._get_span_vec(self._zeta)
        self.u_a_cp_unswept = self._get_unswept_axial_vec(self._zeta)
        self.u_n_cp_unswept = self._get_unswept_normal_vec(self._zeta)
        self.u_s_cp_unswept = self._get_unswept_span_vec(self._zeta)
        self.c_bar_cp = self._get_cp_avg_chord_lengths()
        self.twist_cp = self.get_twist(self._zeta)
        self.dihedral_cp = self.get_dihedral(self._zeta)
        self.sweep_cp = self.get_sweep(self._zeta)
        self.dS = abs(self.node_span_locs[1:]-self.node_span_locs[:-1])*self.b*self.c_bar_cp


    def _airfoil_interpolator(self, interp_spans, sample_spans, coefs):
        # Interpolates the airfoil coefficients at the given span locations.
        # Allows for the coefficients having been evaluated as a function of 
        # span as well.
        # Solution found on stackoverflow
        i = np.arange(interp_spans.size)
        j = np.searchsorted(sample_spans, interp_spans) - 1
        j[np.where(j<0)] = 0 # Not allowed to go outside the array
        d = (interp_spans-sample_spans[j])/(sample_spans[j+1]-sample_spans[j])
        return (1-d)*coefs[i,j]+d*coefs[i,j+1]


    def get_cp_CLa(self, alpha, Rey, Mach):
        """Returns the lift slope at each control point.

        Parameters
        ----------
        alpha : ndarray
            Angle of attack

        Rey : ndarray
            Reynolds number

        Mach : ndarray
            Mach number

        Returns
        -------
        float
            Lift slope
        """

        # Gather lift slopes
        CLas = np.zeros((self._N,self._num_airfoils))
        for j in range(self._num_airfoils):
            CLas[:,j] = self._airfoils[j].get_CLa(alpha=alpha, Rey=Rey, Mach=Mach, trailing_flap_deflection=self._delta_flap, trailing_flap_fraction=self._cp_c_f)

        # Interpolate
        if self._num_airfoils == 1:
            return CLas.flatten()
        else:
            return self._airfoil_interpolator(self._zeta, self._airfoil_spans, CLas)


    def get_cp_aL0(self, Rey, Mach):
        """Returns the zero-lift angle of attack at each control point.

        Parameters
        ----------
        Rey : ndarray
            Reynolds number

        Mach : ndarray
            Mach number

        Returns
        -------
        float
            Zero lift angle of attack
        """

        # Gather zero-lift angles of attack
        aL0s = np.zeros((self._N,self._num_airfoils))
        for j in range(self._num_airfoils):
            aL0s[:,j] = self._airfoils[j].get_aL0(Rey=Rey, Mach=Mach, trailing_flap_deflection=self._delta_flap, trailing_flap_fraction=self._cp_c_f)

        # Interpolate
        if self._num_airfoils == 1:
            return aL0s.flatten()
        else:
            return self._airfoil_interpolator(self._zeta, self._airfoil_spans, aL0s)


    def get_cp_CLRe(self, alpha, Rey, Mach):
        """Returns the derivative of the lift coefficient with respect to Reynolds number at each control point

        Parameters
        ----------
        alpha : ndarray
            Angle of attack

        Rey : ndarray
            Reynolds number

        Mach : ndarray
            Mach number

        Returns
        -------
        float
            Z
        """

        # Gather Reynolds slopes
        CLRes = np.zeros((self._N,self._num_airfoils))
        for j in range(self._num_airfoils):
            CLRes[:,j] = self._airfoils[j].get_CLRe(alpha=alpha, Rey=Rey, Mach=Mach, trailing_flap_deflection=self._delta_flap, trailing_flap_fraction=self._cp_c_f)

        # Interpolate
        if self._num_airfoils == 1:
            return CLRes.flatten()
        else:
            return self._airfoil_interpolator(self._zeta, self._airfoil_spans, CLRes)

    
    def get_cp_CLM(self, alpha, Rey, Mach):
        """Returns the derivative of the lift coefficient with respect to Mach number at each control point

        Parameters
        ----------
        alpha : ndarray
            Angle of attack

        Rey : ndarray
            Reynolds number

        Mach : ndarray
            Mach number

        Returns
        -------
        float
            Z
        """

        # Get Mach slopes
        CLMs = np.zeros((self._N,self._num_airfoils))
        for j in range(self._num_airfoils):
            CLMs[:,j] = self._airfoils[j].get_CLM(alpha=alpha, Rey=Rey, Mach=Mach, trailing_flap_deflection=self._delta_flap, trailing_flap_fraction=self._cp_c_f)

        # Interpolate
        if self._num_airfoils == 1:
            return CLMs.flatten()
        else:
            return self._airfoil_interpolator(self._zeta, self._airfoil_spans, CLMs)


    def get_cp_CL(self, alpha, Rey, Mach):
        """Returns the coefficient of lift at each control point as a function of params.

        Parameters
        ----------
        alpha : ndarray
            Angle of attack

        Rey : ndarray
            Reynolds number

        Mach : ndarray
            Mach number

        Returns
        -------
        float or ndarray
            Coefficient of lift
        """

        # Get CL
        CLs = np.zeros((self._N,self._num_airfoils))
        for j in range(self._num_airfoils):
            CLs[:,j] = self._airfoils[j].get_CL(alpha=alpha, Rey=Rey, Mach=Mach, trailing_flap_deflection=self._delta_flap, trailing_flap_fraction=self._cp_c_f)

        # Interpolate
        if self._num_airfoils == 1:
            return CLs.flatten()
        else:
            return self._airfoil_interpolator(self._zeta, self._airfoil_spans, CLs)


    def get_cp_CD(self, alpha, Rey, Mach):
        """Returns the coefficient of drag at each control point as a function of params.

        Parameters
        ----------
        alpha : ndarray
            Angle of attack

        Rey : ndarray
            Reynolds number

        Mach : ndarray
            Mach number

        Returns
        -------
        float
            Coefficient of drag
        """

        # Get CD
        CDs = np.zeros((self._N,self._num_airfoils))
        for j in range(self._num_airfoils):
            CDs[:,j] = self._airfoils[j].get_CD(alpha=alpha, Rey=Rey, Mach=Mach, trailing_flap_deflection=self._delta_flap, trailing_flap_fraction=self._cp_c_f)

        # Interpolate
        if self._num_airfoils == 1:
            return CDs.flatten()
        else:
            return self._airfoil_interpolator(self._zeta, self._airfoil_spans, CDs)


    def get_cp_Cm(self, alpha, Rey, Mach):
        """Returns the moment coefficient at each control point as a function of params.

        Parameters
        ----------
        alpha : ndarray
            Angle of attack

        Rey : ndarray
            Reynolds number

        Mach : ndarray
            Mach number

        Returns
        -------
        float
            Moment coefficient
        """

        # Get Cm
        Cms = np.zeros((self._N,self._num_airfoils))
        for j in range(self._num_airfoils):
            Cms[:,j] = self._airfoils[j].get_Cm(alpha=alpha, Rey=Rey, Mach=Mach, trailing_flap_deflection=self._delta_flap, trailing_flap_fraction=self._cp_c_f)

        # Interpolate
        if self._num_airfoils == 1:
            return Cms.flatten()
        else:
            return self._airfoil_interpolator(self._zeta, self._airfoil_spans, Cms)


    def get_outline_points(self):
        """Returns a set of points that represents the planar outline of the blade.
        
        Returns
        -------
        ndarray
            Array of outline points.
        """
        zetas = np.linspace(self._zeta_hub, 1.0, 3*self._N)
        chords = self.get_chord(zetas)
        twists = self.get_twist(zetas)

        # Spanwise locations
        points = np.zeros((self._N*6+1,2))
        points[:self._N,0] = zetas*self.diameter
        points[-2:self._N-1:-1,1] = zetas*self.diameter

        # Leading edge
        points[:self._N,1] = 0.25*np.cos(twists)*chords

        # Trailing edge
        points[-2:self._N-1:-1,1] = -0.75*np.cos(twists)*chords

        # Complete the circle
        points[-1,:] = points[0,:]
        return points


    def get_stl_vectors(self, **kwargs):
        """Calculates and returns the outline vectors required for 
        generating an .stl model of the wing segment.

        Parameters
        ----------
        section_resolution : int, optional
            Number of points to use in distcretizing the airfoil sections. Defaults to 200.

        close_te : bool, optional
            Whether to force the trailing edge to be sealed. Defaults to true

        Returns
        -------
        ndarray
            Array of outline vectors. First index is the facet index, second is the point
            index, third is the vector components.
        """

        # Discretize by node locations
        section_res = kwargs.get("section_resolution", 200)
        close_te = kwargs.get("close_te", True)
        num_facets = self._N*(section_res-1)*2
        vectors = np.zeros((num_facets*3,3))

        # Generate vectors
        for i in range(self._N):

            # Root-ward node
            root_span = self.node_span_locs[i]
            root_outline = self._get_airfoil_outline_coords_at_span(root_span, section_res, close_te)

            # Tip-ward node
            tip_span = self.node_span_locs[i+1]
            tip_outline = self._get_airfoil_outline_coords_at_span(tip_span, section_res, close_te)

            # Create facets between the outlines
            for j in range(section_res-1):
                index = (2*i*(section_res-1)+2*j)*3

                vectors[index] = root_outline[j]
                vectors[index+1] = tip_outline[j+1]
                vectors[index+2] = tip_outline[j]

                vectors[index+3] = tip_outline[j+1]
                vectors[index+4] = root_outline[j]
                vectors[index+5] = root_outline[j+1]

        return vectors


    def _get_airfoil_outline_coords_at_span(self, span, _N, close_te):
        # Returns the airfoil section outline in body-fixed coordinates at the specified span fraction with the specified number of points

        # Determine flap deflection and fraction at this point
        if self._has_control_surface and span >= self._cntrl_root_span and span <= self._cntrl_tip_span:
            if self.side == "left":
                d_f = np.interp(span, self._zeta[::-1], self._delta_flap[::-1])
            else:
                d_f = np.interp(span, self._zeta, self._delta_flap)
            c_f = self.get_c_f(span)
        else:
            d_f = 0.0
            c_f = 0.0

        # Linearly interpolate outlines, ignoring twist, etc for now
        if self._num_airfoils == 1:
            points = self._airfoils[0].get_outline_points(_N=_N, trailing_flap_deflection=d_f, trailing_flap_fraction=c_f, close_te=close_te)
        else:
            index = 0
            while True:
                if span >= self._airfoil_spans[index] and span <= self._airfoil_spans[index+1]:
                    total_span = self._airfoil_spans[index+1]-self._airfoil_spans[index]

                    # Get weights
                    root_weight = 1-abs(span-self._airfoil_spans[index])/total_span
                    tip_weight = 1-abs(span-self._airfoil_spans[index+1])/total_span

                    # Get outlines
                    root_outline = self._airfoils[index].get_outline_points(_N=_N, trailing_flap_deflection=d_f, trailing_flap_fraction=c_f, close_te=close_te)
                    tip_outline = self._airfoils[index+1].get_outline_points(_N=_N, trailing_flap_deflection=d_f, trailing_flap_fraction=c_f, close_te=close_te)

                    # Interpolate
                    points = root_weight*root_outline+tip_weight*tip_outline
                    break

                index += 1

        # Get twist, dihedral, and chord
        twist = self.get_twist(span)
        dihedral = self.get_dihedral(span)
        chord = self.get_chord(span)

        # Scale to chord and transform to body-fixed coordinates
        if self._shear_dihedral:
            q = euler_to_quat(np.array([0.0, twist, 0.0]))
        else:
            q = euler_to_quat(np.array([dihedral, twist, 0.0]))
        untransformed_coords = chord*np.array([-points[:,0].flatten()+0.25, np.zeros(_N), -points[:,1]]).T
        coords = self._get_quarter_chord_loc(span)[np.newaxis]+quat_inv_trans(q, untransformed_coords)

        # Seal trailing edge
        te = (coords[0]+coords[-1])*0.5
        coords[0] = te
        coords[-1] = te

        return coords


    def export_stp(self, **kwargs):
        """Creates a FreeCAD part representing a loft of the wing segment.

        Parameters
        ----------
        airplane_name: str
            Name of the airplane this segment belongs to.

        file_tag : str, optional
            Optional tag to prepend to output filename default. The output files will be named "<AIRCRAFT_NAME>_<WING_NAME>.stp".

        section_resolution : int
            Number of outline points to use for the sections. Defaults to 200.
        
        spline : bool, optional
            Whether the wing segment sections should be represented using splines. This can cause issues with some geometries/CAD 
            packages. Defaults to False.

        maintain_sections : bool, optional
            Whether the wing segment sections should be preserved in the loft. Defaults to True.

        close_te : bool, optional
            Whether to force the trailing edge to be sealed. Defaults to true
        """

        # Import necessary modules
        import FreeCAD
        import Part

        # Kwargs
        airplane_name = kwargs.get("airplane_name")
        file_tag = kwargs.get("file_tag", "")
        section_resolution = kwargs.get("section_resolution", 200)
        spline = kwargs.get("spline", False)
        maintain_sections = kwargs.get("maintain_sections", True)
        close_te = kwargs.get("close_te", True)

        # Create sections
        sections = []
        for s_i in self.node_span_locs:
            points = []

            # Get outline points
            outline = self._get_airfoil_outline_coords_at_span(s_i, section_resolution, close_te)

            # Check for wing going to a point
            if np.all(np.all(outline == outline[0,:])):
                #tip = FreeCAD.Base.Vector(*outline[0])
                #points.append(tip)
                #continue
                #TODO loft to an actual point
                outline = self._get_airfoil_outline_coords_at_span(s_i-0.000001, section_resolution, close_te)

            # Create outline points
            for point in outline:
                points.append(FreeCAD.Base.Vector(*point))

            # Add to section list
            if not spline: # Use polygon
                section_polygon = Part.makePolygon(points)
                sections.append(section_polygon)
            else: # Use spline
                section_spline = Part.BSplineCurve(points)
                sections.append(section_spline.toShape())

        # Loft
        wing_loft = Part.makeLoft(sections, True, maintain_sections, False).Faces
        wing_shell = Part.Shell(wing_loft)
        wing_solid = Part.Solid(wing_shell)

        # Export
        abs_path = os.path.abspath("{0}{1}_{2}.stp".format(file_tag, airplane_name, self.name))
        wing_solid.exportStep(abs_path)


    def export_dxf(self, airplane_name, **kwargs):
        """Creates a dxf representing successive sections of the wing segment.

        Parameters
        ----------
        airplane_name: str
            Name of the airplane this segment belongs to.

        file_tag : str, optional
            Optional tag to prepend to output filename default. The output files will be named "<AIRCRAFT_NAME>_<WING_NAME>.stp".

        section_resolution : int
            Number of outline points to use for the sections. Defaults to 200.

        close_te : bool, optional
            Whether to force the trailing edge to be sealed. Defaults to true
        """

        # Get kwargs
        file_tag = kwargs.get("file_tag", "")
        section_res = kwargs.get("section_resolution", 200)
        close_te = kwargs.get("close_te", True)

        # Initialize arrays
        X = np.zeros((self._N+1, section_res))
        Y = np.zeros((self._N+1, section_res))
        Z = np.zeros((self._N+1, section_res))

        # Fill arrays
        for i, s_i in enumerate(self.node_span_locs):

            # Get outline points
            outline = self._get_airfoil_outline_coords_at_span(s_i, section_res, close_te)

            # Store in arrays
            X[i,:] = outline[:,0]
            Y[i,:] = outline[:,1]
            Z[i,:] = outline[:,2]

        # Export
        abs_path = os.path.abspath("{0}{1}_{2}.dxf".format(file_tag, airplane_name, self.name))
        dxf_spline(abs_path, X, Y, Z)