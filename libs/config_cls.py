import ROOT
import os

class Configuration:
    def __init__(self,  baseline, sig_file_path, bkg_file_path, tree_name, 
                        pdf_dir = './pdf', datacard_dir = './datacards'):
        self.baseline      = baseline
        self.sig_file_path = sig_file_path
        self.bkg_file_path = bkg_file_path
        self.tree_name     = tree_name
        self.pdf_dir       = pdf_dir             
        self.categories = []

        if not os.path.exists(pdf_dir)      : os.mkdir(pdf_dir)
        if not os.path.exists(datacard_dir) : os.mkdir(datacard_dir)

    def add_category(self, category):
        print '[INFO] adding category', category.name, '%s' %"BLINDED" if category.blind else "NOT BLINDED"

        category.set_sig_file_path(self.sig_file_path)
        category.set_bkg_file_path(self.bkg_file_path)
        category.set_tree_name    (self.tree_name    )
        category.set_pdf_dir      (self.pdf_dir      )
        category.update_selection (self.baseline     )

        self.categories.append(category)
    
    def fit_model(self):
        for cc in self.categories: cc.fit_model()
        