import io
import math

import numpy as np
import pandas as pd
import rdkit
import rdkit.Chem
import rdkit.Chem.Draw
import streamlit as st


# Sample data is a subset of DRD2 dataset.
sample_data = [
    "CCCC=N[SH](=O)(O)C=Cc1ccccc1",
    "O=C(Nc1ccc2c(c1)OCO2)c1cnn2cccnc12",
    "N#Cc1c(N)n(CCNC(=O)c2cccs2)c2nc3ccccc3nc12",
    "CC(=O)c1cccc(NC(=O)Cn2cnc3c2c(=O)n(C)c(=O)n3C)c1",
    "CCOC(=O)c1oc2ccccc2c1NC(=O)Cc1cccs1",
    "CN1CCCC(OC(=O)C(C#N)=Cc2ccccc2)C1",
    "CC(C)N(Cc1ccccc1)C(=S)NC1CCCCC1",
    "Cc1ccccc1OCCC(=O)OCC(=O)NC1CCCC1",
    "COc1ccc(C(CNC(=O)c2cccc(NS(=O)(=O)c3ccc(C)c(F)c3)c2)N2CCOCC2)cc1",
    "CC(C)(C)N(NC(=O)Nc1ccc(Cl)c(Cl)c1)C(=O)c1ccccc1",
    "CC(=O)Nc1ccc(OCC(O)CN2C(C)(C)CC(O)CC2(C)C)cc1",
    "COc1cccc(NC(=O)C(Sc2nnnn2C)c2ccccc2)c1",
    "O=C(O)CCCCCn1c(SC=C2CC(=O)N3C=C(Cl)C=CC3=N2)nc2ccsc2c1=O",
    "COc1ccc(N(CC(=O)NCc2ccccc2OC)S(=O)(=O)c2ccc(C)cc2)cc1Cl",
]


def rdimage(mol: rdkit.Chem.Mol) -> bytes:
    """Returns an image of a molecule as bytes."""
    pil_image = rdkit.Chem.Draw.MolToImage(mol)
    buf = io.BytesIO()
    pil_image.save(buf, format="png")
    return buf.getvalue()


@st.cache
def rdimage_from_smi(smi: str) -> bytes:
    return rdimage(rdkit.Chem.MolFromSmiles(smi))


def find_row_idx(i):
    """Returns 'i'-th element of st.session_state.df that has 'rating' of 0 (Unscored).

    Non-zero rating indicates that molecule was rated already,
    and we would like to omit rated molecules from display.
    """

    df = st.session_state.df
    zs = np.where(df["rating"].to_numpy() == 0)[0]
    if i < len(zs):
        return zs[i]
    else:
        return None


def on_change(j):
    """Callback for radio widget to record selected rating into df."""
    st.session_state.df.at[j, "rating"] = st.session_state[f"rating_{j}"]


def make_grid():

    df = st.session_state.df

    N_COLS = st.sidebar.number_input(label="Number of columns", value=3, min_value=1)
    N_ROWS = math.ceil(len(df) / N_COLS)  # Rows to fit all molecules.

    for row in range(N_ROWS):

        # We want to display radio widget side-by-side with the molecule.
        # Without columns Streamlit displays elements under each other.
        # Streamlit columns could be used, but we use them for the grid,
        # and Streamlit does not allow to nest columns.
        # To get side-by-side effect, create twice as many columns,
        # and use even {0, 2...} for mols and odd {1, 3...} for radio widget.
        cols2 = st.columns(N_COLS * 2)

        for col in range(N_COLS):

            # We want to omit molecules that were rated already.
            # For that, we use find_row_idx.
            # Set j = i for displaying all molecules, even rated ones.
            i = row * N_COLS + col  # GUI "Cell" index.
            j = find_row_idx(i)  # Pandas dataframe index.

            if j is not None:
                left_col = cols2[col * 2]  # For mol image.
                right_col = cols2[col * 2 + 1]  # For radio widget.
                with left_col:
                    st.image(
                        rdimage_from_smi(df["SMILES"][j]),
                        # caption=df["SMILES"][j],
                        use_column_width="always",
                    )
                with right_col:
                    radio_choices = ["Unscored"] + list(range(1, 6))
                    value = st.radio(
                        label="",
                        options=range(6),  # Internal representation, 0,1,2,3,4,5.
                        format_func=lambda c: radio_choices[c],  # User-visible labels.
                        key=f"rating_{j}",
                        on_change=on_change,
                        args=(j,),
                    )


st.set_page_config(
    page_title="MolWall",
    layout="wide",
)


# Set to None for interactive upload, set to str to read from a local file.
local_file_name = None


def add_df_to_session(df):
    st.session_state["df"] = df
    st.session_state["df"]["rating"] = 0


# We store dataframe in streamlit.session_state.
# This way we don't read it on every redraw,
# and we can store additional data in the dataframe (column "rating").
if "df" not in st.session_state:
    if local_file_name is not None:
        df = pd.read_csv(local_file_name)
        add_df_to_session(df)
    else:
        st.write("""
            ### Welcome to Mol Wall.
            
            Please upload a CSV file with a SMILES column 
            (see [examples](https://github.com/Augmented-Drug-Design-Human-in-the-Loop/mol-wall/blob/main/tests/data/)),
            or load sample molecules.
        """)
        label = "Upload a CSV file with SMILES column"
        uploaded_file = st.sidebar.file_uploader(label, type=["csv", "csv.gz"])
        if uploaded_file is not None:
            df = pd.read_csv(uploaded_file)
            add_df_to_session(df)
            st.session_state["uploaded_file"] = uploaded_file.name

        btn_result = st.sidebar.button("...or load sample molecules")
        if btn_result:
            df = pd.DataFrame({"SMILES": sample_data})
            add_df_to_session(df)
            st.session_state["uploaded_file"] = None

if "df" in st.session_state:
    make_grid()

    st.download_button(
        label="Download data as CSV",
        data=st.session_state.df.to_csv(index=False).encode("utf-8"),
        file_name=st.session_state["uploaded_file"].split('.csv')[0] + "_with_ratings.csv"  if st.session_state["uploaded_file"] is not None else "file_with_ratings.csv",
        mime="text/csv",
    )
