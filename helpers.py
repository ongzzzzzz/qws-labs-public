import copy
import numpy as np
import itertools

# Bokeh interactive plotting
from bokeh.plotting import figure, output_notebook, show
from bokeh.models import Slider, ColumnDataSource, CustomJS, Range1d, Toggle, Button, SetValue
from bokeh.layouts import column, row
from bokeh.io import push_notebook

def convert_to_basis_string(arr):

    arr = copy.deepcopy(arr)
    arr = np.array(arr, dtype=object)
    arr[arr == 0] = 'Z'
    arr[arr == 1] = 'X'

    return arr

def get_all_pauli_strings(length):
	return [''.join(p) for p in list(itertools.product("IXYZ", repeat=length))]


def interactive_sine_waves():

    # Create the x-values for the sine wave (ranging from 0 to 2*pi)
    x = np.linspace(0, 2 * np.pi, 500)
    
    # Initialize the data source for sine waves
    source = ColumnDataSource(data=dict(x=x, y1=np.sin(x), y2=np.sin(x), y3=2*np.sin(x)))
    
    # Create a plot
    plot = figure(title="Interactive Sine Waves", x_axis_label='x', y_axis_label='y', width=800, height=400)
    
    # Add the sine wave lines to the plot
    plot.line('x', 'y1', source=source, line_width=2, line_dash = [3], color="blue", legend_label="Sine Wave 1")
    plot.line('x', 'y2', source=source, line_width=2, line_dash = [5], color="red", legend_label="Sine Wave 2")
    plot.line('x', 'y3', source=source, line_width=2, color="green", legend_label="Sine Wave 1 + Sine Wave 2")
    
    plot.y_range = Range1d(-5, 5)
    
    # Create sliders for amplitude, frequency, and phase for each sine wave
    amplitude_slider1 = Slider(start=0, end=2.3, value=1, step=0.1, title="Amplitude #1 [a.u.]")
    frequency_slider1 = Slider(start=0.1, end=2, value=1, step=0.1, title="Frequency #1 [1/s]")
    phase_slider1 = Slider(start=0, end=360, value=0, step=5, title="Phase #1 [deg.]")
    
    amplitude_slider2 = Slider(start=0, end=2.3, value=1, step=0.1, title="Amplitude #2 [a.u.]")
    frequency_slider2 = Slider(start=0.1, end=2, value=1, step=0.1, title="Frequency #2 [1/s]")
    phase_slider2 = Slider(start=0, end=360, value=0, step=5, title="Phase #2 [deg.]")
    
    # CustomJS callback to update sine waves and their sum
    callback = CustomJS(args=dict(source=source, x=x, amplitude_slider1=amplitude_slider1,
                                  frequency_slider1=frequency_slider1, phase_slider1=phase_slider1,
                                  amplitude_slider2=amplitude_slider2, frequency_slider2=frequency_slider2,
                                  phase_slider2=phase_slider2), code="""
        var data = source.data;
        var amplitude1 = amplitude_slider1.value;
        var frequency1 = frequency_slider1.value;
        var phase1 = phase_slider1.value;
        
        var amplitude2 = amplitude_slider2.value;
        var frequency2 = frequency_slider2.value;
        var phase2 = phase_slider2.value;
        
        var y1 = [];
        var y2 = [];
        var y3 = [];
        
        for (var i = 0; i < x.length; i++) {
            y1.push(amplitude1 * Math.sin(frequency1 * x[i] + (Math.PI/180)*phase1));
            y2.push(amplitude2 * Math.sin(frequency2 * x[i] + (Math.PI/180)*phase2));
            y3.push(y1[i] + y2[i]);  // Sum of the two sine waves
        }
        
        data['y1'] = y1;
        data['y2'] = y2;
        data['y3'] = y3;  // Update the third sine wave (sum)
        source.change.emit();
    """)
    
    # Attach the CustomJS callback to the sliders
    amplitude_slider1.js_on_change('value', callback)
    frequency_slider1.js_on_change('value', callback)
    phase_slider1.js_on_change('value', callback)
    
    amplitude_slider2.js_on_change('value', callback)
    frequency_slider2.js_on_change('value', callback)
    phase_slider2.js_on_change('value', callback)
    
    # Layout the plot and sliders in a column
    layout = column(row( column(phase_slider1, amplitude_slider1, frequency_slider1), column(phase_slider2, amplitude_slider2, frequency_slider2) ), plot)
    
    # Show the plot in the notebook
    show(layout)

class double_slit_experiment():
    
    def __init__(self, d=10e-6, l=2, wavelength=50e-12):
        # Initialization of parameters
        self.d = d
        self.l = l
        self.wavelength = wavelength
        
        self.num_particles = 1000
        
        # Set the x axis and related values for simulation
        self.set_x_axis()
        
        # Histogram
        bins = np.linspace(self.x[0], self.x[-1], 50)
        hist, edges = np.histogram([], density=False, bins=bins)
        
#         self.source = ColumnDataSource(data=dict(y=hist, x=hist, left=, right=))
        
        # Initialize the figure and axis for plotting
        self.hist_source = ColumnDataSource(data=dict(start=edges[:-1], end=edges[1:], top=hist))  # Data for histogram
        self.geometry_source = ColumnDataSource(data=dict(d=[self.d], l=[self.l], wavelength=[self.wavelength], num_particles=[self.num_particles]))

        # Create histogram plot
        self.fig = figure(title="Double Slit Simulation", 
                          x_axis_label='Distance from center [m]',
                          y_axis_label='Number of Electrons',
                          width=800, height=400)

        self.fig.x_range = Range1d(-5e-5, 5e-5)
        
        self.fig.quad(top='top', bottom=0, left='start', right='end', source=self.hist_source, 
                      fill_color="blue", line_color="white", alpha=0.6)
        
        # Create a Button widget
        self.toggle_button = Toggle(label="Ignore", button_type="primary", width=100)

        tb_callback = CustomJS(args=dict(tb=self.toggle_button), code="""
            const a = cb_obj.label;
             if (a === 'Observe') {
                tb.label = "Ignore";
            } else {
                tb.label = 'Observe';
            }
        """)

        # JavaScript callback to toggle the label between '0' and '1'
        self.toggle_button.js_on_click(tb_callback)
        
        # Create a Button widget to fire electrons
        self.button = Button(label="Fire Electrons", button_type="success")
        
        # Create CustomJS callback to fire particles
        self.button.js_on_click(self.create_custom_callback())
        
        callback_d = CustomJS(args=dict(source=self.geometry_source), code="""
            const f = cb_obj.value 
            source.data['d'][0] = f * 1e-6
        """)

        self.d_slider = Slider(start=5, end=20, value=self.d*1e6, step=.1, title=r"Distance between slits (um)")
        self.d_slider.js_on_change('value', callback_d)
        
        callback_l = CustomJS(args=dict(source=self.geometry_source), code="""
            const f = cb_obj.value 
            source.data['l'][0] = f
        """)

        self.l_slider = Slider(start=1, end=3, value=self.l, step=.1, title=r"Distance to screen (m)")
        self.l_slider.js_on_change('value', callback_l)
        
        callback_num_particles = CustomJS(args=dict(source=self.geometry_source), code="""
            const f = cb_obj.value 
            source.data['num_particles'][0] = f
        """)

        self.num_particles_slider = Slider(start=5, end=2000, value=self.num_particles, step=5, title=r"Number of electrons")
        self.num_particles_slider.js_on_change('value', callback_num_particles)

        # Show the layout in the notebook
        show(self.layout(), notebook_handle=True)
    
        
    def set_x_axis(self, width_scale=5, points=1000):
        # Generate x-axis points and corresponding angles
        self.x = np.linspace(-width_scale * self.d, width_scale * self.d, points)
        theta = np.arctan(self.x / self.l)
        self.delta_l = self.d * np.sin(theta)
    
    def create_custom_callback(self):
        # CustomJS callback that fires particles and updates the plot
        return CustomJS(args=dict(hist_source=self.hist_source, x=self.x, geometry_source=self.geometry_source, tb=self.toggle_button), code="""
            var num_particles = geometry_source.data['num_particles'][0];
            var pmf = [];
            
            const d = geometry_source.data['d'][0];
            const l = geometry_source.data['l'][0];
            const wavelength = geometry_source.data['wavelength'][0];
            
            var delta_l = [];
            for (var i = 0; i < x.length; i++) {
                delta_l.push(d * Math.sin(x[i] / l));
            }
            
            const observe_status = tb.label
            
            // Calculate the pdf
            var fringes = (function() {
                var intensity = [];
                var gaussian = [];
                for (var i = 0; i < x.length; i++) {
                    var cos_value = Math.cos(Math.PI * delta_l[i] / wavelength);
                    if (observe_status === 'Ignore') {
                        intensity.push(cos_value * cos_value);
                    } else {
                        intensity.push(1);
                    }
                    gaussian.push(Math.exp(-0.5 * Math.pow(x[i] / (2 * d), 2)) / Math.exp(0));
                }

                var fringes = [];
                for (var i = 0; i < x.length; i++) {
                    fringes.push(gaussian[i] * intensity[i]);
                }

                // Normalize
                var sum = fringes.reduce(function(a, b) { return a + b; }, 0);
                for (var i = 0; i < fringes.length; i++) {
                    fringes[i] /= sum;
                }

                return fringes;
            })();
            
            // Sample particles from the probability mass function
            var result = [];
            for (var i = 0; i < num_particles; i++) {
                var rand = Math.random();
                var cumulative_prob = 0;
                for (var j = 0; j < x.length; j++) {
                    cumulative_prob += fringes[j];
                    if (rand < cumulative_prob) {
                        result.push(x[j]);
                        break;
                    }
                }
            }
           
            // Calculate the histogram data based on the new particle positions
            var hist_bins = 200;  // Number of histogram bins
            var hist_start = Math.min.apply(null, result);
            var hist_end = Math.max.apply(null, result);
            var hist_step = (hist_end - hist_start) / hist_bins;
            var bins = new Array(hist_bins).fill(0);
            
            // Compute the histogram
            for (var i = 0; i < result.length; i++) {
                var index = Math.floor((result[i] - hist_start) / hist_step);
                if (index >= 0 && index < hist_bins) {
                    bins[index]++;
                }
            }
            
            // Prepare the new histogram data
            var start_vals = [];
            var end_vals = [];
            var top_vals = [];
            for (var i = 0; i < hist_bins; i++) {
                start_vals.push((hist_start + i * hist_step)/1);
                end_vals.push((hist_start + (i + 1) * hist_step)/1);
                top_vals.push(bins[i]);
            }

            // Update the histogram data source
            hist_source.data['start'] = start_vals;
            hist_source.data['end'] = end_vals;
            hist_source.data['top'] = top_vals;
            hist_source.change.emit();  // Notify Bokeh to update the histogram
        """)
    
    def layout(self):
        # Return the layout that includes the button and plot
        return column(self.d_slider,self.l_slider, self.num_particles_slider, row(self.toggle_button,self.button), self.fig)
