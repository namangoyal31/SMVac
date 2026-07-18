import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def create_contour_plot():
    print("Loading data...")
    df = pd.read_csv('results/numerical_data.csv')

    stable = df[df['Stability'] == 1]
    metastable = df[df['Stability'] == 2]
    unstable = df[df['Stability'] == 3]
    nonpert = df[df['Stability'] == 4]

    plt.figure(figsize=(10, 8))
    
    # Plot the background scatter points
    plt.scatter(stable['Mh_calc'], stable['Mt'], c='green', alpha=0.3, s=5, label='Stable')
    plt.scatter(metastable['Mh_calc'], metastable['Mt'], c='yellow', alpha=0.3, s=5, label='Metastable')
    plt.scatter(unstable['Mh_calc'], unstable['Mt'], c='red', alpha=0.3, s=5, label='Unstable')
    
    if not nonpert.empty:
        plt.scatter(nonpert['Mh_calc'], nonpert['Mt'], c='black', alpha=0.3, s=5, label='Non-perturbative')

    # To prevent contours from bleeding into the stable region, we MUST NOT delete the stable points.
    # But using all 1 million stable points makes Delaunay triangulation too slow.
    # To prevent contours from bleeding into the stable region, we MUST NOT delete the stable points.
    # But using all 1 million stable points can be slow. Let's use all active points and 10% of stable points.
    df_active = df[df['Stability'] != 1].copy()
    df_stable = df[df['Stability'] == 1].sample(frac=0.1, random_state=42).copy()
    
    df_contour = pd.concat([df_active, df_stable], ignore_index=True)
    
    # Handle any Infs or NaNs
    df_contour.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_contour.dropna(subset=['S_exact', 'Mh_calc', 'Mt'], inplace=True)
    
    df_contour.loc[df_contour['S_exact'] < 0, 'S_exact'] = 1e8
    
    # Clamp classical instability artifacts (huge positive values) to 0, as they are deeply unstable
    df_contour.loc[df_contour['S_exact'] > 1e10, 'S_exact'] = 0.0
    
    if not df_contour.empty:
        print("Drawing contours...")
        
        # Define contour levels to span the Unstable (S < 450) and Metastable (S > 450) regions
        # We add 500000 specifically to draw a contour through the Mh=5, Mt=96 point.
        levels = [10, 50, 100, 200, 450, 1000, 3000, 500000]
        
        # tricontour handles irregularly spaced data (due to Mh_calc being calculated dynamically)
        contour = plt.tricontour(df_contour['Mh_calc'], df_contour['Mt'], df_contour['S_exact'], 
                                 levels=levels, colors='blue', linewidths=1.5, alpha=0.9)
        
        # Add labels to the contour lines, manually placing them along Mh=50
        # We explicitly define the y-coordinates to perfectly hit each of the 8 contours exactly once
        y_locs = [127, 140, 155, 170, 185, 200, 215, 230]
        manual_locations = [(50, y) for y in y_locs]
        
        plt.clabel(contour, inline=True, fontsize=11, fmt='S=%1.0f', manual=manual_locations)

    plt.scatter(125.10, 173.1, color='blue', marker='*', s=200, edgecolor='black', label='Experimental (125.1, 173.1)', zorder=5)

    plt.xlabel('Higgs Mass $M_h$ (GeV)', fontsize=12)
    plt.ylabel('Top Mass $M_t$ (GeV)', fontsize=12)
    plt.title('Numerical Phase Diagram with Constant Action Contours', fontsize=14)
    
    # Set the same bounds as the standard numerical plot (0 to 250)
    plt.xlim(0, 250)
    plt.ylim(0, 250)
        
    plt.legend(loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    output_file = 'results/numerical_stability_contours.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {output_file}")

if __name__ == '__main__':
    create_contour_plot()
