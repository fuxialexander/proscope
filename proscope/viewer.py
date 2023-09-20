# %%
# use nglview to visualize the structure
import nglview as nv
import py3Dmol

from atac_rna_data_processing.io.s3_utils import open_file_with_s3

def view_pdb_ipython(pdb_path):
    view = nv.NGLWidget()
    view.add_component(nv.FileStructure(pdb_path))
    # view.update_cartoon(color='red', component=0)
    return view

def view_pdb_html(path_to_pdb, s3_file_sys=None):
    '''
        #function to display pdb in py3dmol
    '''
    pdb = open_file_with_s3(path_to_pdb, "r", s3_file_sys=s3_file_sys).read()
    
    view = py3Dmol.view(width=500, height=500)
    view.addModel(pdb, "pdb")
    # color by chain

    view.setStyle({'chain': 'A'}, {'cartoon': {'color': 'red'}})
    view.setStyle({'chain': 'B'}, {'cartoon': {'color': 'blue'}})
    # Add hydrogen bonds between chains A and B (or any specific atoms)
    # view.addHBonds({'chain': 'A'}, {'chain': 'B'}, {'dist': 3.5, 'color': 'yellow', 'linewidth': 3})

    view.zoomTo()
    output = view._make_html().replace("'", '"')
    x = f"""<!DOCTYPE html><html></center> {output} </center></html>"""  # do not use ' in this input
    
    return f"""<iframe height="500px" width="100%"  name="result" allow="midi; geolocation; microphone; camera;
                            display-capture; encrypted-media;" sandbox="allow-modals allow-forms
                            allow-scripts allow-same-origin allow-popups
                            allow-top-navigation-by-user-activation allow-downloads" allowfullscreen=""
                            allowpaymentrequest="" frameborder="0" srcdoc='{x}'></iframe>"""

'''
    return f"""<iframe  style="width: 100%; height:700px" name="result" allow="midi; geolocation; microphone; camera; 
                            display-capture; encrypted-media;" sandbox="allow-modals allow-forms 
                            allow-scripts allow-same-origin allow-popups 
                            allow-top-navigation-by-user-activation allow-downloads" allowfullscreen="" 
                            allowpaymentrequest="" frameborder="0" srcdoc='{x}'></iframe>"""
'''


