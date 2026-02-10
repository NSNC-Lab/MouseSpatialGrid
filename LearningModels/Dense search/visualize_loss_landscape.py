"""
5D Loss Landscape Visualization Suite
=====================================
Visualizes the loss landscape from a 5D grid search.

Data format:
- losses: (epochs, score_type, grid_points) - shape (64, 2, 16384)
- params: (epochs, parameters, grid_points) - shape (64, 5, 16384)

Grid structure: The 16384 points come from varying 4 parameters 
(one param is constant), creating a 4×16×16×16 = 16384 grid.
"""

import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from itertools import combinations
import os

# Try to import optional dependencies
try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Warning: Plotly not installed. Some interactive plots will be unavailable.")
    print("Install with: pip install plotly")

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False
    print("Warning: Pandas not installed. Some plots will be unavailable.")
    print("Install with: pip install pandas")


class LossLandscapeVisualizer:
    """
    Visualizer for 5D loss landscape from grid search results.
    """
    
    def __init__(self, mat_file_path, score_type=1):
        """
        Initialize the visualizer.
        
        Parameters:
        -----------
        mat_file_path : str
            Path to the .mat file containing 'losses' and 'params'
        score_type : int
            0 for L2, 1 for MSE (default: 1 for MSE)
        """
        self.mat_file_path = mat_file_path
        self.score_type = score_type
        self.score_name = "L2" if score_type == 0 else "MSE"
        
        # Load data
        self._load_data()
        
        # Setup output directory
        self.output_dir = os.path.join(os.path.dirname(mat_file_path), "LossLandscapePlots")
        os.makedirs(self.output_dir, exist_ok=True)
        
    def _load_data(self):
        """Load and process the .mat file data."""
        data = sio.loadmat(self.mat_file_path)
        
        # The data has shape (batches, ..., points_per_batch)
        # Flatten across batches to get the full grid search
        raw_params = data['params']  # (batches, 5, points_per_batch)
        raw_losses = data['losses']  # (batches, 2, points_per_batch)
        
        n_batches = raw_params.shape[0]
        n_points_per_batch = raw_params.shape[2]
        total_points = n_batches * n_points_per_batch
        
        print(f"Data structure: {n_batches} batches × {n_points_per_batch} points = {total_points:,} total grid points")
        
        # Flatten: concatenate all batches into one continuous array
        # params: (5, total_points), losses: (total_points,)
        self.params = raw_params.transpose(1, 0, 2).reshape(5, -1)  # (5, total_points)
        self.losses = raw_losses[:, self.score_type, :].flatten()  # (total_points,)
        
        self.best_params = data.get('best_params', None)
        self.best_loss = data.get('best_loss', None)
        
        # Extract unique parameter values and names
        self.param_names = ['Param 0', 'Param 1', 'Param 2', 'Param 3', 'Param 4']
        self.param_values = []
        self.param_sizes = []
        
        for i in range(5):
            unique_vals = np.unique(self.params[i])
            self.param_values.append(unique_vals)
            self.param_sizes.append(len(unique_vals))
            
        print(f"Loaded {len(self.losses):,} grid points")
        print(f"Parameter grid sizes: {self.param_sizes}")
        print(f"Expected from grid: {np.prod(self.param_sizes):,}")
        print(f"Score type: {self.score_name}")
        
        # Reshape flat array into grid if possible
        self._reshape_to_grid()
        
        # Find optimal point
        self.optimal_flat_idx = np.argmin(self.losses)
        self.optimal_loss = self.losses[self.optimal_flat_idx]
        self.optimal_params = self.params[:, self.optimal_flat_idx]
        
        print(f"\nOptimal {self.score_name} loss: {self.optimal_loss:.4f}")
        print(f"Optimal params: {self.optimal_params}")
        
    def _reshape_to_grid(self):
        """
        Attempt to reshape flat grid points into a multidimensional grid.
        This requires understanding how the grid was flattened.
        """
        # Try to infer the grid structure from the parameter ordering
        n_points = len(self.losses)
        
        # Check if we can reshape based on unique value counts
        expected_size = np.prod(self.param_sizes)
        
        if expected_size == n_points:
            # Build mapping from parameter values to grid indices
            self.grid_shape = tuple(self.param_sizes)
            self.loss_grid = np.full(self.grid_shape, np.nan)
            self.param_to_idx = []
            
            for i in range(5):
                mapping = {v: idx for idx, v in enumerate(self.param_values[i])}
                self.param_to_idx.append(mapping)
            
            # Fill the grid
            for flat_idx in range(n_points):
                grid_idx = tuple(
                    self.param_to_idx[i][self.params[i, flat_idx]] 
                    for i in range(5)
                )
                self.loss_grid[grid_idx] = self.losses[flat_idx]
            
            # Find optimal grid index
            self.optimal_grid_idx = np.unravel_index(np.nanargmin(self.loss_grid), self.grid_shape)
            print(f"Reshaped to grid: {self.grid_shape}")
            print(f"Optimal grid index: {self.optimal_grid_idx}")
            self.has_grid = True
        else:
            print(f"Warning: Cannot reshape to regular grid. Expected {expected_size}, got {n_points}")
            self.has_grid = False
            self.loss_grid = None
            
    def set_param_names(self, names):
        """Set custom parameter names."""
        if len(names) != 5:
            raise ValueError("Must provide exactly 5 parameter names")
        self.param_names = names
        
    # =========================================================================
    # VISUALIZATION 1: 1D Slices
    # =========================================================================
    def plot_1d_slices(self, save=True, show=True, eprop_params=None):
        """
        Plot 1D slices through the loss landscape.
        Fix 4 parameters at optimal values, vary one parameter at a time.
        
        Parameters:
        -----------
        eprop_params : list or array, optional
            Reference parameter set to highlight (e.g., E-prop optimal)
        """
        if not self.has_grid:
            print("Cannot create 1D slices without regular grid structure.")
            return None
        
        # Default E-prop optimal parameters
        if eprop_params is None:
            eprop_params = [0.0062, 5.36e-04, 0.0081, 0, 0.1355]
            
        # Only plot params that have more than 1 unique value
        active_params = [i for i in range(5) if self.param_sizes[i] > 1]
        n_active = len(active_params)
        
        fig, axes = plt.subplots(1, n_active, figsize=(4*n_active, 4))
        if n_active == 1:
            axes = [axes]
        
        for ax_idx, i in enumerate(active_params):
            # Create slicer that fixes all params at optimal except param i
            slicer = list(self.optimal_grid_idx)
            slicer[i] = slice(None)
            slice_1d = self.loss_grid[tuple(slicer)]
            
            ax = axes[ax_idx]
            ax.plot(self.param_values[i], slice_1d, 'b-o', linewidth=2, markersize=6)
            ax.axvline(self.optimal_params[i], color='r', linestyle='--', 
                      linewidth=2, label='Optimal Grid')
            ax.axvline(eprop_params[i], color='blue', linestyle=':', 
                      linewidth=2, label='Optimal E-prop')
            ax.axhline(self.optimal_loss, color='g', linestyle=':', 
                      alpha=0.5, label='Min Loss')
            ax.set_xlabel(self.param_names[i], fontsize=12)
            ax.set_ylabel(f'{self.score_name} Loss', fontsize=12)
            ax.set_title(f'Loss vs {self.param_names[i]}', fontsize=12)
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)
        
        plt.suptitle(f'1D Slices Through {self.score_name} Loss Landscape', fontsize=14)
        plt.tight_layout()
        
        if save:
            path = os.path.join(self.output_dir, 'loss_landscape_1d_slices.png')
            plt.savefig(path, dpi=150, bbox_inches='tight')
            print(f"Saved: {path}")
            
        if show:
            plt.show()
        else:
            plt.close()
            
        return fig
    
    # =========================================================================
    # VISUALIZATION 2: 2D Pairwise Heatmaps
    # =========================================================================
    def plot_2d_slices(self, save=True, show=True, log_scale=False):
        """
        Plot 2D heatmap slices for all pairs of parameters.
        Fix 3 parameters at optimal, vary 2 at a time.
        """
        if not self.has_grid:
            print("Cannot create 2D slices without regular grid structure.")
            return None
        
        # Only consider params with more than 1 unique value
        active_params = [i for i in range(5) if self.param_sizes[i] > 1]
        pairs = list(combinations(active_params, 2))
        n_pairs = len(pairs)
        
        if n_pairs == 0:
            print("Not enough varying parameters for 2D plots.")
            return None
        
        # Determine subplot layout
        n_cols = min(3, n_pairs)
        n_rows = (n_pairs + n_cols - 1) // n_cols
        
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
        if n_pairs == 1:
            axes = np.array([[axes]])
        axes = np.atleast_2d(axes)
        
        for idx, (i, j) in enumerate(pairs):
            row, col = idx // n_cols, idx % n_cols
            ax = axes[row, col]
            
            # Create slicer that fixes all params except i and j
            slicer = list(self.optimal_grid_idx)
            slicer[i] = slice(None)
            slicer[j] = slice(None)
            slice_2d = self.loss_grid[tuple(slicer)]
            
            # The slice_2d has shape (size_i, size_j)
            # For imshow, rows=y-axis (param i), cols=x-axis (param j)
            
            plot_data = np.log10(slice_2d) if log_scale else slice_2d
            
            im = ax.imshow(plot_data, origin='lower', aspect='auto',
                          extent=[self.param_values[j][0], self.param_values[j][-1],
                                  self.param_values[i][0], self.param_values[i][-1]],
                          cmap='viridis')
            
            # Mark optimal point
            ax.scatter([self.optimal_params[j]], [self.optimal_params[i]], 
                      c='red', s=100, marker='*', edgecolors='white', linewidths=1,
                      label='Optimal', zorder=5)
            
            ax.set_xlabel(self.param_names[j], fontsize=11)
            ax.set_ylabel(self.param_names[i], fontsize=11)
            ax.set_title(f'{self.param_names[i]} vs {self.param_names[j]}', fontsize=11)
            
            cbar_label = f'log10({self.score_name} Loss)' if log_scale else f'{self.score_name} Loss'
            cbar = plt.colorbar(im, ax=ax, label=cbar_label)
        
        # Hide unused subplots
        for idx in range(n_pairs, n_rows * n_cols):
            row, col = idx // n_cols, idx % n_cols
            axes[row, col].axis('off')
        
        scale_str = " (Log Scale)" if log_scale else ""
        plt.suptitle(f'2D {self.score_name} Loss Landscape Slices{scale_str}', fontsize=14)
        plt.tight_layout()
        
        if save:
            suffix = "_log" if log_scale else ""
            path = os.path.join(self.output_dir, f'loss_landscape_2d_slices{suffix}.png')
            plt.savefig(path, dpi=150, bbox_inches='tight')
            print(f"Saved: {path}")
            
        if show:
            plt.show()
        else:
            plt.close()
            
        return fig
    
    # =========================================================================
    # VISUALIZATION 3: Parallel Coordinates Plot
    # =========================================================================
    def plot_parallel_coordinates(self, save=True, show=True, n_samples=None):
        """
        Create a parallel coordinates plot colored by loss.
        Shows all 5 parameters + loss as parallel axes.
        
        Parameters:
        -----------
        n_samples : int, optional
            Number of samples to plot (for performance). If None, plot all.
        """
        if not PLOTLY_AVAILABLE or not PANDAS_AVAILABLE:
            print("Parallel coordinates requires plotly and pandas.")
            return None
        
        # Create DataFrame
        df = pd.DataFrame({
            name: self.params[i] for i, name in enumerate(self.param_names)
        })
        df['Loss'] = self.losses
        
        # Subsample if requested
        if n_samples is not None and n_samples < len(df):
            df = df.sample(n=n_samples, random_state=42)
        
        # Only include params that vary
        active_dims = [name for i, name in enumerate(self.param_names) 
                       if self.param_sizes[i] > 1]
        
        fig = px.parallel_coordinates(
            df, 
            dimensions=active_dims + ['Loss'],
            color='Loss',
            color_continuous_scale='Viridis_r',  # Reversed so lower loss = brighter
            title=f'Parallel Coordinates Plot - {self.score_name} Loss ({len(df):,} points)'
        )
        
        fig.update_layout(
            font=dict(size=12),
            coloraxis_colorbar=dict(title=f'{self.score_name} Loss')
        )
        
        if save:
            path = os.path.join(self.output_dir, 'loss_landscape_parallel_coords.html')
            fig.write_html(path)
            print(f"Saved: {path}")
            
        if show:
            fig.show()
            
        return fig
    
    # =========================================================================
    # VISUALIZATION 4: Interactive 2D Slice with Sliders
    # =========================================================================
    def plot_interactive_2d_slice(self, fixed_params=(0, 1, 2), vary_params=(3, 4)):
        """
        Create an interactive 2D heatmap with sliders for the fixed parameters.
        
        Parameters:
        -----------
        fixed_params : tuple
            Indices of 3 parameters to control with sliders
        vary_params : tuple
            Indices of 2 parameters to show in the heatmap
        """
        if not self.has_grid:
            print("Cannot create interactive plot without regular grid structure.")
            return None
        
        # Validate inputs
        all_params = set(fixed_params) | set(vary_params)
        if len(all_params) != 5:
            print("Must specify all 5 parameters between fixed_params and vary_params")
            return None
            
        fig, ax = plt.subplots(figsize=(10, 8))
        plt.subplots_adjust(bottom=0.35)
        
        i, j = vary_params  # y-axis, x-axis
        
        # Initial indices for fixed params (middle values)
        current_idx = {p: self.param_sizes[p] // 2 for p in fixed_params}
        
        def get_slice():
            slicer = [None] * 5
            slicer[i] = slice(None)
            slicer[j] = slice(None)
            for p in fixed_params:
                slicer[p] = current_idx[p]
            return self.loss_grid[tuple(slicer)]
        
        # Initial plot
        slice_2d = get_slice()
        im = ax.imshow(slice_2d, origin='lower', aspect='auto',
                      extent=[self.param_values[j][0], self.param_values[j][-1],
                              self.param_values[i][0], self.param_values[i][-1]],
                      cmap='viridis')
        
        ax.set_xlabel(self.param_names[j], fontsize=12)
        ax.set_ylabel(self.param_names[i], fontsize=12)
        ax.set_title(f'{self.param_names[i]} vs {self.param_names[j]}')
        cbar = plt.colorbar(im, ax=ax, label=f'{self.score_name} Loss')
        
        # Create sliders for fixed parameters
        sliders = {}
        slider_axes = {}
        
        for slider_idx, p in enumerate(fixed_params):
            if self.param_sizes[p] <= 1:
                continue  # Skip constant parameters
                
            ax_slider = plt.axes([0.2, 0.25 - slider_idx * 0.07, 0.6, 0.04])
            slider = Slider(
                ax_slider, 
                self.param_names[p],
                0, 
                self.param_sizes[p] - 1,
                valinit=current_idx[p],
                valstep=1
            )
            sliders[p] = slider
            slider_axes[p] = ax_slider
            
            # Add value label
            ax_slider.text(1.05, 0.5, f'{self.param_values[p][current_idx[p]]:.4f}',
                          transform=ax_slider.transAxes, va='center')
        
        def update(val):
            for p, slider in sliders.items():
                current_idx[p] = int(slider.val)
                # Update value label
                slider.ax.texts[-1].set_text(f'{self.param_values[p][current_idx[p]]:.4f}')
            
            new_slice = get_slice()
            im.set_data(new_slice)
            im.set_clim(new_slice.min(), new_slice.max())
            fig.canvas.draw_idle()
        
        for slider in sliders.values():
            slider.on_changed(update)
        
        plt.suptitle(f'Interactive {self.score_name} Loss Landscape', fontsize=14, y=0.98)
        plt.show()
        
        return fig
    
    # =========================================================================
    # VISUALIZATION 5: Marginal Distribution Plots
    # =========================================================================
    def plot_marginal_losses(self, save=True, show=True):
        """
        Plot marginal loss distributions.
        Shows mean loss (marginalized over other params) with min/max range.
        """
        if not self.has_grid:
            print("Cannot create marginal plots without regular grid structure.")
            return None
        
        active_params = [i for i in range(5) if self.param_sizes[i] > 1]
        n_active = len(active_params)
        
        fig, axes = plt.subplots(1, n_active, figsize=(4*n_active, 4))
        if n_active == 1:
            axes = [axes]
        
        for ax_idx, i in enumerate(active_params):
            ax = axes[ax_idx]
            
            # Compute marginal statistics over all other dimensions
            axes_to_reduce = tuple(j for j in range(5) if j != i)
            
            marginal_mean = np.nanmean(self.loss_grid, axis=axes_to_reduce)
            marginal_min = np.nanmin(self.loss_grid, axis=axes_to_reduce)
            marginal_max = np.nanmax(self.loss_grid, axis=axes_to_reduce)
            marginal_std = np.nanstd(self.loss_grid, axis=axes_to_reduce)
            
            x = self.param_values[i]
            
            # Plot mean with std band
            ax.plot(x, marginal_mean, 'b-o', linewidth=2, markersize=6, label='Mean')
            ax.fill_between(x, marginal_mean - marginal_std, marginal_mean + marginal_std,
                           alpha=0.3, color='blue', label='±1 Std')
            
            # Plot min/max range
            ax.fill_between(x, marginal_min, marginal_max, alpha=0.1, color='gray',
                           label='Min-Max Range')
            ax.plot(x, marginal_min, 'g--', alpha=0.5, linewidth=1)
            ax.plot(x, marginal_max, 'r--', alpha=0.5, linewidth=1)
            
            # Mark optimal
            ax.axvline(self.optimal_params[i], color='red', linestyle=':', 
                      linewidth=2, label='Optimal')
            
            ax.set_xlabel(self.param_names[i], fontsize=12)
            ax.set_ylabel(f'{self.score_name} Loss', fontsize=12)
            ax.set_title(f'Marginal: {self.param_names[i]}', fontsize=12)
            ax.legend(loc='upper right', fontsize=8)
            ax.grid(True, alpha=0.3)
        
        plt.suptitle(f'Marginal {self.score_name} Loss Distributions', fontsize=14)
        plt.tight_layout()
        
        if save:
            path = os.path.join(self.output_dir, 'loss_landscape_marginals.png')
            plt.savefig(path, dpi=150, bbox_inches='tight')
            print(f"Saved: {path}")
            
        if show:
            plt.show()
        else:
            plt.close()
            
        return fig
    
    # =========================================================================
    # VISUALIZATION 6: Loss Distribution
    # =========================================================================
    def plot_loss_distribution(self, save=True, show=True):
        """
        Plot the distribution of losses across all grid points.
        """
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))
        
        # Left: Histogram
        ax = axes[0]
        ax.hist(self.losses, bins=100, edgecolor='black', alpha=0.7, color='steelblue')
        ax.axvline(self.optimal_loss, color='red', linestyle='--', 
                  linewidth=2, label=f'Best: {self.optimal_loss:.2f}')
        ax.axvline(np.median(self.losses), color='orange', linestyle=':', 
                  linewidth=2, label=f'Median: {np.median(self.losses):.2f}')
        ax.set_xlabel(f'{self.score_name} Loss', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title(f'Distribution of {self.score_name} Losses ({len(self.losses):,} points)', fontsize=12)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        # Right: Sorted losses (like a loss landscape profile)
        ax = axes[1]
        sorted_losses = np.sort(self.losses)
        percentiles = np.linspace(0, 100, len(sorted_losses))
        ax.plot(percentiles, sorted_losses, 'b-', linewidth=1)
        ax.axhline(self.optimal_loss, color='red', linestyle='--', 
                  linewidth=2, label=f'Best: {self.optimal_loss:.2f}')
        
        # Mark key percentiles
        for p in [1, 5, 10, 25, 50]:
            val = np.percentile(self.losses, p)
            ax.axhline(val, color='gray', linestyle=':', alpha=0.5)
            ax.text(p + 2, val, f'{p}%: {val:.1f}', fontsize=9, va='bottom')
        
        ax.set_xlabel('Percentile', fontsize=12)
        ax.set_ylabel(f'{self.score_name} Loss', fontsize=12)
        ax.set_title('Loss Percentile Profile', fontsize=12)
        ax.legend(loc='upper left')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 100)
        
        plt.suptitle(f'{self.score_name} Loss Distribution Analysis', fontsize=14)
        plt.tight_layout()
        
        if save:
            path = os.path.join(self.output_dir, 'loss_distribution.png')
            plt.savefig(path, dpi=150, bbox_inches='tight')
            print(f"Saved: {path}")
            
        if show:
            plt.show()
        else:
            plt.close()
            
        return fig
    
    # =========================================================================
    # VISUALIZATION 7: 3D Surface Plot (Plotly)
    # =========================================================================
    def plot_3d_surface(self, param_x=3, param_y=4, save=True, show=True):
        """
        Create an interactive 3D surface plot for two parameters.
        
        Parameters:
        -----------
        param_x, param_y : int
            Indices of parameters for x and y axes
        """
        if not PLOTLY_AVAILABLE:
            print("3D surface plot requires plotly.")
            return None
            
        if not self.has_grid:
            print("Cannot create 3D surface without regular grid structure.")
            return None
        
        # Get 2D slice at optimal values for other params
        slicer = list(self.optimal_grid_idx)
        slicer[param_x] = slice(None)
        slicer[param_y] = slice(None)
        slice_2d = self.loss_grid[tuple(slicer)]
        
        # Create meshgrid for plotting
        X, Y = np.meshgrid(self.param_values[param_y], self.param_values[param_x])
        
        fig = go.Figure(data=[go.Surface(
            x=X, y=Y, z=slice_2d,
            colorscale='Viridis',
            colorbar=dict(title=f'{self.score_name} Loss')
        )])
        
        # Add optimal point marker
        opt_x = self.optimal_params[param_y]
        opt_y = self.optimal_params[param_x]
        fig.add_trace(go.Scatter3d(
            x=[opt_x], y=[opt_y], z=[self.optimal_loss],
            mode='markers',
            marker=dict(size=10, color='red', symbol='diamond'),
            name='Optimal'
        ))
        
        fig.update_layout(
            title=f'3D {self.score_name} Loss Landscape',
            scene=dict(
                xaxis_title=self.param_names[param_y],
                yaxis_title=self.param_names[param_x],
                zaxis_title=f'{self.score_name} Loss'
            ),
            width=900,
            height=700
        )
        
        if save:
            path = os.path.join(self.output_dir, 'loss_landscape_3d_surface.html')
            fig.write_html(path)
            print(f"Saved: {path}")
            
        if show:
            fig.show()
            
        return fig
    
    # =========================================================================
    # VISUALIZATION 8: Contour Plot with Gradient Arrows
    # =========================================================================
    def plot_contour_with_gradient(self, param_x=3, param_y=4, save=True, show=True):
        """
        Create a contour plot with gradient direction arrows.
        """
        if not self.has_grid:
            print("Cannot create contour plot without regular grid structure.")
            return None
        
        # Get 2D slice
        slicer = list(self.optimal_grid_idx)
        slicer[param_x] = slice(None)
        slicer[param_y] = slice(None)
        slice_2d = self.loss_grid[tuple(slicer)]
        
        # Create meshgrid
        X, Y = np.meshgrid(self.param_values[param_y], self.param_values[param_x])
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Contour plot
        levels = 20
        contour = ax.contourf(X, Y, slice_2d, levels=levels, cmap='viridis')
        ax.contour(X, Y, slice_2d, levels=levels, colors='white', alpha=0.3, linewidths=0.5)
        plt.colorbar(contour, ax=ax, label=f'{self.score_name} Loss')
        
        # Compute and plot gradients
        grad_y, grad_x = np.gradient(slice_2d)
        
        # Subsample for cleaner arrows
        step = max(1, min(slice_2d.shape) // 8)
        ax.quiver(X[::step, ::step], Y[::step, ::step], 
                 -grad_x[::step, ::step], -grad_y[::step, ::step],
                 color='white', alpha=0.7, scale=None)
        
        # Mark optimal
        ax.scatter([self.optimal_params[param_y]], [self.optimal_params[param_x]], 
                  c='red', s=200, marker='*', edgecolors='white', linewidths=2,
                  label='Optimal', zorder=5)
        
        ax.set_xlabel(self.param_names[param_y], fontsize=12)
        ax.set_ylabel(self.param_names[param_x], fontsize=12)
        ax.set_title(f'{self.score_name} Loss Contours with Gradient Directions', fontsize=14)
        ax.legend()
        
        plt.tight_layout()
        
        if save:
            path = os.path.join(self.output_dir, 'loss_landscape_contour_gradient.png')
            plt.savefig(path, dpi=150, bbox_inches='tight')
            print(f"Saved: {path}")
            
        if show:
            plt.show()
        else:
            plt.close()
            
        return fig
    
    # =========================================================================
    # VISUALIZATION 9: All 3D Surface Combinations
    # =========================================================================
    def generate_all_3d_surfaces(self, save=True, show=False, combined=True):
        """
        Generate 3D surface plots for all combinations of 2 parameters.
        With 5 parameters, this creates C(5,2) = 10 plots.
        
        Parameters:
        -----------
        combined : bool
            If True, save all surfaces in a single HTML file.
            If False, save each surface as a separate HTML file.
        """
        if not PLOTLY_AVAILABLE:
            print("3D surface plots require plotly.")
            return None
            
        if not self.has_grid:
            print("Cannot create 3D surfaces without regular grid structure.")
            return None
        
        # Create subfolder for 3D surfaces
        surfaces_dir = os.path.join(self.output_dir, "3D_Surfaces")
        os.makedirs(surfaces_dir, exist_ok=True)
        
        # Get all pairs of parameters (that have more than 1 unique value)
        active_params = [i for i in range(5) if self.param_sizes[i] > 1]
        pairs = list(combinations(active_params, 2))
        n_pairs = len(pairs)
        
        print(f"Generating {n_pairs} 3D surface plots...")
        
        if combined:
            # Create a single HTML with all surfaces using subplots
            n_cols = 2
            n_rows = (n_pairs + n_cols - 1) // n_cols
            
            # Create subplot specs for 3D scenes
            specs = [[{'type': 'scene'} for _ in range(n_cols)] for _ in range(n_rows)]
            
            subplot_titles = []
            for param_x, param_y in pairs:
                subplot_titles.append(f'{self.param_names[param_x]} vs {self.param_names[param_y]}')
            
            fig = make_subplots(
                rows=n_rows, cols=n_cols,
                specs=specs,
                subplot_titles=subplot_titles,
                horizontal_spacing=0.05,
                vertical_spacing=0.08
            )
            
            for idx, (param_x, param_y) in enumerate(pairs):
                row = idx // n_cols + 1
                col = idx % n_cols + 1
                scene_name = f'scene{idx + 1}' if idx > 0 else 'scene'
                
                # Get 2D slice at optimal values for other params
                slicer = list(self.optimal_grid_idx)
                slicer[param_x] = slice(None)
                slicer[param_y] = slice(None)
                slice_2d = self.loss_grid[tuple(slicer)]
                
                # Create meshgrid for plotting
                X, Y = np.meshgrid(self.param_values[param_y], self.param_values[param_x])
                
                # Add surface
                fig.add_trace(
                    go.Surface(
                        x=X, y=Y, z=slice_2d,
                        colorscale='Viridis',
                        showscale=(idx == 0),  # Only show colorbar for first plot
                        name=f'{self.param_names[param_x]} vs {self.param_names[param_y]}'
                    ),
                    row=row, col=col
                )
                
                # Add optimal point marker
                opt_x = self.optimal_params[param_y]
                opt_y = self.optimal_params[param_x]
                opt_z = self.optimal_loss
                fig.add_trace(
                    go.Scatter3d(
                        x=[opt_x], y=[opt_y], z=[opt_z],
                        mode='markers',
                        marker=dict(size=5, color='red', symbol='diamond'),
                        showlegend=False
                    ),
                    row=row, col=col
                )
                
                # Update scene axes
                fig.update_scenes(
                    dict(
                        xaxis_title=self.param_names[param_y],
                        yaxis_title=self.param_names[param_x],
                        zaxis_title=self.score_name
                    ),
                    row=row, col=col
                )
                
                print(f"  [{idx+1}/{n_pairs}] Added: {self.param_names[param_x]} vs {self.param_names[param_y]}")
            
            fig.update_layout(
                title=f'{self.score_name} Loss Landscape - All Parameter Combinations',
                height=400 * n_rows,
                width=1400,
                showlegend=False
            )
            
            if save:
                path = os.path.join(surfaces_dir, 'all_3d_surfaces_combined.html')
                fig.write_html(path)
                print(f"Saved combined file: {path}")
            
            if show:
                fig.show()
                
            return fig
        
        else:
            # Save each surface as a separate file
            figures = []
            for idx, (param_x, param_y) in enumerate(pairs):
                # Get 2D slice at optimal values for other params
                slicer = list(self.optimal_grid_idx)
                slicer[param_x] = slice(None)
                slicer[param_y] = slice(None)
                slice_2d = self.loss_grid[tuple(slicer)]
                
                # Create meshgrid for plotting
                X, Y = np.meshgrid(self.param_values[param_y], self.param_values[param_x])
                
                fig = go.Figure(data=[go.Surface(
                    x=X, y=Y, z=slice_2d,
                    colorscale='Viridis',
                    colorbar=dict(title=f'{self.score_name} Loss')
                )])
                
                # Add optimal point marker
                opt_x = self.optimal_params[param_y]
                opt_y = self.optimal_params[param_x]
                opt_z = self.optimal_loss
                fig.add_trace(go.Scatter3d(
                    x=[opt_x], y=[opt_y], z=[opt_z],
                    mode='markers',
                    marker=dict(size=8, color='red', symbol='diamond'),
                    name='Optimal'
                ))
                
                name_x = self.param_names[param_x]
                name_y = self.param_names[param_y]
                
                fig.update_layout(
                    title=f'{self.score_name} Loss: {name_x} vs {name_y}',
                    scene=dict(
                        xaxis_title=name_y,
                        yaxis_title=name_x,
                        zaxis_title=f'{self.score_name} Loss'
                    ),
                    width=800,
                    height=600
                )
                
                if save:
                    # Create safe filename
                    safe_name_x = name_x.replace(' ', '_').replace('/', '_')
                    safe_name_y = name_y.replace(' ', '_').replace('/', '_')
                    filename = f'3d_surface_{safe_name_x}_vs_{safe_name_y}.html'
                    path = os.path.join(surfaces_dir, filename)
                    fig.write_html(path)
                    print(f"  [{idx+1}/{n_pairs}] Saved: {filename}")
                
                if show:
                    fig.show()
                    
                figures.append(fig)
            
            print(f"All 3D surfaces saved to: {surfaces_dir}")
            return figures
    
    # =========================================================================
    # RUN ALL VISUALIZATIONS
    # =========================================================================
    def generate_all_plots(self, show=False):
        """Generate all static visualization plots."""
        print("\n" + "="*60)
        print("Generating all loss landscape visualizations...")
        print("="*60 + "\n")
        
        print("1. 1D Slices...")
        self.plot_1d_slices(save=True, show=show)
        
        print("\n2. 2D Pairwise Heatmaps...")
        self.plot_2d_slices(save=True, show=show)
        self.plot_2d_slices(save=True, show=show, log_scale=True)
        
        print("\n3. Marginal Distributions...")
        self.plot_marginal_losses(save=True, show=show)
        
        print("\n4. Loss Distribution...")
        self.plot_loss_distribution(save=True, show=show)
        
        print("\n5. Contour with Gradients...")
        self.plot_contour_with_gradient(save=True, show=show)
        
        if PLOTLY_AVAILABLE and PANDAS_AVAILABLE:
            print("\n6. Parallel Coordinates (Interactive HTML)...")
            self.plot_parallel_coordinates(save=True, show=show)
            
            print("\n7. 3D Surface (Interactive HTML)...")
            self.plot_3d_surface(save=True, show=show)
            
            print("\n8. All 3D Surface Combinations...")
            self.generate_all_3d_surfaces(save=True, show=show, combined=True)
        
        print("\n" + "="*60)
        print(f"All plots saved to: {self.output_dir}")
        print("="*60)


# =============================================================================
# MAIN EXECUTION
# =============================================================================
if __name__ == "__main__":
    # Path to your .mat file
    mat_file = r"c:\Users\ipboy\Documents\GitHub\MouseSpatialGrid\LearningModels\Dense search\Results\output_compressed_DenseSearch_1.1milPointsAroundEPropSolution.mat"
    
    # Create visualizer (score_type=1 for MSE, score_type=0 for L2)
    viz = LossLandscapeVisualizer(mat_file, score_type=1)
    
    # Optional: Set custom parameter names
    viz.set_param_names(['on_input->output', 'off_input->output', 'on_input->PV', 'off_input->PV', 'PV->output'])
    
    # Generate all static plots (set show=True to display interactively)
    viz.generate_all_plots(show=True)
    
    # Uncomment to run interactive 2D slice with sliders:
    viz.plot_interactive_2d_slice(fixed_params=(0, 1, 2), vary_params=(3, 4))
