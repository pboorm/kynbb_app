"""
Streamlit app to show AGN spectrum as function of logNH
Copyright: Peter Boorman (2021)
"""

import plotly.express as px
import streamlit as st
from scipy.interpolate import RegularGridInterpolator
import pandas as pd
import tqdm
import numpy as np
from astropy.io import fits
from matplotlib.pyplot import get_cmap
from matplotlib.colors import Normalize, rgb2hex
from itertools import product

def enumerated_product(*args):
    """
    https://stackoverflow.com/questions/56430745/enumerating-a-tuple-of-indices-with-itertools-product
    """
    yield from zip(product(*(range(len(x)) for x in args)), product(*args))

def load_table_model(tab_path):
    """
    Function to load an Xspec table model with RegularGridInterpolator
    """
    hdul = fits.open(tab_path)
    data = hdul[1].data

    par_values = {}
    for parname, Nvals, parvals in zip(data["NAME"], data["NUMBVALS"], data["VALUE"]):
        if parname == "nH":
            par_values[parname] = np.log10(parvals[0:Nvals])+22.
        else:
            par_values[parname] = parvals[0:Nvals]
            
        print(parname, par_values[parname])#.min(), par_values[parname].max())

    par_spec = {}    
    Ngridpoints = np.prod([len(pars) for name, pars in par_values.items()])
    print("%s loaded, %d grid points with param names: %s" %(tab_path, Ngridpoints, "; ".join(list(par_values.keys()))))

    energies = hdul[2].data
    E_keV = np.array([0.5*(x1+x2) for (x1, x2) in energies])
    E_keV_widths = np.array([x2-x1 for (x1, x2) in energies])
    
    flux_shape = [len(pars) for parname, pars in par_values.items()]
    flux_shape.append(len(E_keV))
    flux_table = np.empty(shape=flux_shape)
    print("Creating flux array with shape:", flux_shape)

    table_par_values = [val[0] for val in hdul[3].data]
    table_flux_values = [val[1] for val in hdul[3].data]

    with tqdm.tqdm(total = Ngridpoints, position = 0, desc = "Gridding (%s)" %(", ".join(list(par_values.keys())))) as pbar:
        for i, (idx, varpar_values) in enumerate(enumerated_product(*par_values.values())):

            # if (table_par_values[i] == varpar_values).all():
            flux_table[idx] = table_flux_values[i]
            # else:
            #     print(table_par_values[i], varpar_values[i])

            pbar.update(1)

    flux_grid = RegularGridInterpolator(tuple(par_values.values()), flux_table)
    
    return E_keV, E_keV_widths, flux_grid

@st.cache
def get_spec():
    """
    Function to extract the initial spectra
    """

    E_keV, E_keV_widths, kynbb_mod = load_table_model("./kynbb.fits")

    return E_keV, E_keV_widths, kynbb_mod

DEGREE_SYMBOL = "\N{DEGREE SIGN}"
E_keV, E_keV_widths, kynbb_mod = get_spec()

# if st.sidebar.checkbox("E * E * model", False):
factor = E_keV*E_keV
ylabel = "Flux / arb."#"keV<sup>2</sup> (Photons cm<sup>-2</sup> s<sup>-1</sup> keV<sup>-1</sup>)"
yrange = [0., 7.]
ytickvals = [1., 10., 100., 1000., 1.e4, 1.e5, 1.e6, 1.e7]
# else:
#     factor = 1.
#     ylabel = "Photons cm<sup>-2</sup> s<sup>-1</sup> keV<sup>-1</sup>"
#     yrange = [-5., 0.]
#     ytickvals = [1.e-5, 1.e-4, 1.e-3, 1.e-2, 1.e-1, 1.]

spin_c = st.sidebar.slider(
        "Spin",
        min_value=0.,
        max_value=0.998,
        value = 0.9,
        step=0.01,
        format="%.2f",
        key="spin",
    )
incl_c = st.sidebar.slider(
        "Inclination",
        min_value=45.,
        max_value=85.,
        value = 60.,
        step=1.,
        format=f"%.0f{DEGREE_SYMBOL}",
        key="incl",
    )
logMBH_c = st.sidebar.slider(
        "log MBH / Msun",
        min_value=4.,
        max_value=9.,
        value = 6.,
        step=0.1,
        format="%.1f",
        key="logMBH",
    )
logarate_c = st.sidebar.slider(
        "log arate / MEdd",
        min_value=-3.,
        max_value=0.,
        value = -2.,
        step=0.1,
        format="%.1f",
        key="logarate",
    )


df = pd.DataFrame(data = {"E_keV": E_keV,
                          "kynbb": factor*kynbb_mod([spin_c, incl_c, logMBH_c, logarate_c])[0]/E_keV_widths,
                          }
                  )

# df = df.sort_values(by=["E_keV"]).reset_index(drop=True)
# # df.loc[:, "Flux"] = df["kynbb"].sum(axis=1)
# df = df.rename(columns={"kynbb": "Flux"})
# order = ["Flux"]
# df = df[["E_keV"]+order]

mod_df = df.melt(
        id_vars=["E_keV"], value_vars=["kynbb"]
    )

mod_df = mod_df.rename(columns={"variable": "Modcomp", "value": "Flux"})




print(mod_df)


## Construct our plot
fig = px.line(
    mod_df,
    x="E_keV",
    y="Flux",
    log_x=True,
    log_y=True,
    width=800,
    height=600,
    labels=dict(Flux=ylabel, E_keV="Energy / keV"),
    color_discrete_sequence=["#fc6d3c"]
)



fig.update_xaxes(ticks="inside", tickwidth = 2.5, ticklen = 10., linewidth = 2.5, linecolor = "black", mirror = True, gridcolor = "LightGray")
fig.update_yaxes(ticks="inside", tickwidth = 2.5, ticklen = 10., linewidth = 2.5, linecolor = "black", mirror = True, gridcolor = "LightGray")

fig.update_traces(line=dict(
                  width=5.))#{"Transmitted": 3., "Reprocessed": 3., "Scattered": 3., "Total": 6.}))

fig.update_layout(
                  # sliders=sliders,
                  showlegend=False,
                  plot_bgcolor = "rgba(0, 0, 0, 0)",
                  legend=dict(yanchor="top",
                              y=0.99,
                              xanchor="left",
                              x=0.01),
                  yaxis=dict(
                             range=yrange,
                             title_text=ylabel,
                             tickfont = dict(size=20),
                             tickvals=ytickvals,
                             tickmode="array",
                             mirror="allticks",
                             side="top",
                             titlefont=dict(size=30)), 
                  xaxis=dict(
                             range=[np.log10(0.01), np.log10(10.)],
                             title_text="Energy / keV",
                             tickfont = dict(size=20),
                             tickvals=[0.01, 0.1, 1., 10.],
                             tickmode="array",
                             mirror="allticks",
                             side="bottom",
                             titlefont=dict(size=30)))

fig['data'][0]['line']['width']=6.
# fig['data'][1]['line']['width']=5.
# fig['data'][2]['line']['width']=5.
# fig['data'][3]['line']['width']=5.

st.plotly_chart(fig)


# # from astropy.cosmology import FlatLambdaCDM
# # Boorman_cosmo = FlatLambdaCDM(H0=67.3, Om0=0.315)








# fig, ax = plt.subplots(figsize=(12, 9))

# logNHlos = 24.
# PhoIndex = 2.
# Ecut = 100.
# TORsigma = 28.
# CTKcover = 0.
# Theta_inc = 60.
# fscat = 1.e-2
# z = 0.
# norm = 1.


# units = "eem"

# T = uxc_trans([logNHlos, PhoIndex, Ecut, TORsigma, CTKcover, Theta_inc])[0]/E_keV_widths
# R = uxc_refl([logNHlos, PhoIndex, Ecut, TORsigma, CTKcover, Theta_inc])[0]/E_keV_widths
# S = uxc_omni([logNHlos, PhoIndex, Ecut, TORsigma, CTKcover, Theta_inc])[0]/E_keV_widths

# T *= norm
# R *= norm
# S *= norm

# if units == "eem":
#     T *= E_keV**2
#     R *= E_keV**2
#     S *= E_keV**2
#     ylabel = r"keV$^{2}$\,(ph\,cm$^{-2}$\,s$^{-1}$\,keV$^{-1}$)"
# else:
#     ylabel = r"ph\,cm$^{-2}$\,s$^{-1}$\,keV$^{-1}$"
    
# LDcm = Boorman_cosmo.luminosity_distance(z).value*3.846e18*1.e6
# # T = T/(1.+z)#((4*np.pi*LDcm**2)*(1+z))
# # R = R/(1.+z)#((4*np.pi*LDcm**2)*(1+z))
# # S = S/(1.+z)#((4*np.pi*LDcm**2)*(1+z))

# S *= fscat

# tot = T + R + S
# ax.plot(E_keV, tot, c="0.", lw=8., zorder=-100., drawstyle="steps-pre", label=r"Total")
# ax.plot(E_keV, T, c="deepskyblue", lw=3., zorder=-5., drawstyle="steps-pre", label=r"Transmitted")
# ax.plot(E_keV, R, c="tomato", lw=3., zorder=-10., drawstyle="steps-pre", label=r"Reprocessed")
# ax.plot(E_keV, S, c="mediumseagreen", lw=3., zorder=-10., drawstyle="steps-pre", label=r"Scattered")

# ax.set_xscale("log")
# ax.set_yscale("log")

# ax.set_xlim([0.3, 80.])
# ax.set_ylim([np.max(np.hstack((T, R, S)).ravel())/1.e5, np.max(np.hstack((T, R, S)).ravel())* 2.])
# print(np.max(np.hstack((T, R, S)).ravel())/1.e5, np.max(np.hstack((T, R, S)).ravel())* 2.)

# ax.set_xticklabels([r"%.5g" %(t) for t in ax.get_xticks()])

# ax.set_xlabel(r"Energy\,/\,keV $\longrightarrow$", x=0., y=0., ha="left", va="top", fontsize=45.)
# ax.set_ylabel(ylabel, fontsize=45.)# $\longrightarrow$", x=0., y=0., ha="left", va="top")

# ax.legend(loc="lower center", ncol=2, bbox_to_anchor=(0.5, 1.))


# df = pd.read_csv("streamlit_uxclumpy.csv")
# x_range = [np.log10(0.3), np.log10(250.)]
# y_range = [np.log10(2.e-4), np.log10(0.99)]


# ## controller
# # logNH_c = st.slider(
# #         "",
# #         min_value=21.,
# #         max_value=25.95,
# #         value = 24.,
# #         step=0.05,
# #         format="%.2f",
# #         key="",
# #     )

# cols = [c for c in df.columns if "E_keV" not in c]

# # fig = px.line(
# #     df,
# #     x=df["E_keV"],
# #     y=df["lognh_%.2f" %(logNH_c)],
# #     log_x=True,
# #     log_y=True,
# #     width=1000,
# #     height=500,
# #     labels=dict(Flux="EFE / keV s-1 cm-2", Energy="Energy / keV"))#,
#     # animation_frame="lognh_%.2f" %(logNH_c))

# test_df = (df.set_index("E_keV")
#                   .unstack()
#                   .reset_index()
#                   .rename(columns={"level_0": "log N(H)", 0: "EFE"}))

# # print(test_df.shape)
# # print(test_df)
# test_df["log N(H)"] = test_df["log N(H)"].map(lambda x: float(x.split("_")[-1]))

# # print(test_df.shape)
# # print(test_df)
# # colour = "rgba(%s)" %(",".join(["%.5f" %(f) for f in cmap_cols[df.columns.get_loc("lognh_%.2f" %(logNH_c)) - 1]]))
# y_range_lin = [10 ** v for v in y_range]
# x_range_lin = [10 ** v for v in x_range]

# cmap = get_cmap('plasma_r')
# norm = Normalize(vmin = 20., vmax = 26.)
# cmap_cols = cmap(norm(test_df["log N(H)"]))
# color_discrete_map = {lognh: rgb2hex(cmap_cols[i]) for i, lognh in enumerate(test_df["log N(H)"].values)}

# ## for line with colourscale: https://community.plotly.com/t/plotly-express-line-chart-color/27333/4
# fig = px.line(test_df, x="E_keV",
#               y="EFE",
#               animation_frame="log N(H)",
#               color="log N(H)",
#               color_discrete_map=color_discrete_map,
#               log_x = True,
#               log_y=True,
#               range_x=x_range_lin,
#               range_y=y_range_lin,
#               width=1000,
#               height=600,
#               )

# st.subheader("${\\tt UXCLUMPY}$ log $N_{\\rm H}$")

# ## more info: https://plotly.com/python/axes/
# fig.update_xaxes(ticks="inside", tickwidth = 2.5, ticklen = 10., linewidth = 2.5, linecolor = "black", mirror = True, gridcolor = "LightGray")
# fig.update_yaxes(ticks="inside", tickwidth = 2.5, ticklen = 10., linewidth = 2.5, linecolor = "black", mirror = True, gridcolor = "LightGray")

# fig.update_traces(line=dict(
#                   # color=colour,
#                   width=4.))

# fig.update_layout(
#                   # sliders=sliders,
#                   showlegend=False,
#                   plot_bgcolor = "rgba(0, 0, 0, 0)",
#                   legend=dict(yanchor="top",
#                               y=0.99,
#                               xanchor="left",
#                               x=0.01),
#                   yaxis=dict(
#                              # range=y_range,
#                              title_text="keV<sup>2</sup> (ph cm<sup>-2</sup> s<sup>-1</sup> keV<sup>-1</sup>)",
#                              tickfont = dict(size=20),
#                              tickvals=[1.e-3, 1.e-2, 1.e-1],
#                              tickmode="array",
#                              mirror="allticks",
#                              side="top",
#                              titlefont=dict(size=30)), 
#                   xaxis=dict(
#                              # range=x_range,
#                              title_text="Energy / keV",
#                              tickfont = dict(size=20),
#                              tickvals=[1., 10., 100.],
#                              tickmode="array",
#                              mirror="allticks",
#                              side="top",
#                              titlefont=dict(size=30)))

# fig.layout.updatemenus[0].buttons[0].args[1]['frame']['duration'] = 30
# fig.layout.updatemenus[0].buttons[0].args[1]['transition']['duration'] = 5

# ## more here: https://plotly.com/python/configuration-options/
# config = {'staticPlot': True}
# st.plotly_chart(fig, use_container_width=True, config=config)
    

# # st.sidebar.markdown("### Model outputs")
# if st.checkbox("Show Table", False):
#     st.subheader("Raw Data Table")
#     # st.write(df[["E_keV", "lognh_%.2f" %(logNH_c)]], index=False)

# # Some advertising
# st.markdown("[UXCLUMPY](https://github.com/JohannesBuchner/xars/blob/master/doc/uxclumpy.rst) [(Buchner et al., 2019)](https://ui.adsabs.harvard.edu/abs/2019A%26A...629A..16B/abstract), &copy; [Dr. Peter Boorman](https://www.peterboorman.com) & [Dr. Adam Hill](https://www.adambenhill.com)")
