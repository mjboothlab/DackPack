

# ---------------------------------------------------------------------------
# MAIN SCRIPT
# ---------------------------------------------------------------------------
import re
import sys
import PySide6

import pandas as pd

import matplotlib
matplotlib.use("QtAgg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import matplotlib.cm as cm

from PySide6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QLineEdit, QComboBox, QSpinBox, 
    QDoubleSpinBox, QRadioButton, QButtonGroup)

from nucdraw import NucDraw #10.5281/zenodo.15352138
from nupack import Strand, Tube, SetSpec, Model, tube_analysis


# ------------------------------
# STRAND INPUT ROW
# ------------------------------
class StrandRow(QWidget):
    def __init__(self, parent_app):
        super().__init__()
        self.parent_app = parent_app

        layout = QHBoxLayout()
        layout.setSpacing(2)
        layout.setContentsMargins(0, 0, 0, 0)

        self.name = QLineEdit()
        self.name.setPlaceholderText("Name")
        self.name.setFixedWidth(70)

        self.sequence = QLineEdit()
        self.sequence.setPlaceholderText("e.g. ATGCGT")
        self.sequence.setMinimumWidth(140)

        self.conc = QDoubleSpinBox()
        self.conc.setValue(1.00)
        self.conc.setSuffix(" µM")
        self.conc.setFixedWidth(80)

        self.remove_btn = QPushButton("Remove")
        self.remove_btn.setProperty("danger", True)
        self.remove_btn.clicked.connect(self.remove_self)

        layout.addWidget(self.name)
        layout.addWidget(self.sequence)
        layout.addWidget(self.conc)
        layout.addWidget(self.remove_btn)

        self.setLayout(layout)

        self.sequence.textChanged.connect(self.validate_sequence)

    # check only valid characters are entered
    def validate_sequence(self):
        seq = self.sequence.text().upper()
        self.sequence.blockSignals(True)
        self.sequence.setText(seq)
        self.sequence.blockSignals(False)

        valid_bases = "ATGC" if self.parent_app.is_dna() else "AUGC"

        if re.fullmatch(f"[{valid_bases}]*", seq):
            self.sequence.setStyleSheet("")
        else:
            self.sequence.setStyleSheet("border: 2px solid #e84118;")

    # automatically convert between DNA and RNA sequences
    def convert_sequence(self, to_material):
        seq = self.sequence.text().upper()
        if to_material == "RNA":
            seq = seq.replace("T", "U")
        else:
            seq = seq.replace("U", "T")
        self.sequence.setText(seq)
        self.validate_sequence()

    # remove strand row button
    def remove_self(self):
        self.setParent(None)
        self.parent_app.strand_rows.remove(self)


# ------------------------------
# MAIN APP
# ------------------------------
class NupackApp(QWidget):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Chelsea's NUPACK Desktop Tool")
        self.resize(480, 620)

        self.tube_result_global = None
        self.complex_map = {}

        self.setStyleSheet(self.get_styles())

        main_layout = QVBoxLayout()
        main_layout.setSpacing(2)
        main_layout.setContentsMargins(12, 12, 12, 12)
        self.setLayout(main_layout)

        # Section 1
        main_layout.addWidget(self.section_label("1. Model Parameters"))
        main_layout.addWidget(self.model_card())

        # Section 2
        main_layout.addWidget(self.section_label("2. Strand Sequences"))
        main_layout.addWidget(self.strand_card())

        # Section 3
        main_layout.addWidget(self.section_label("3. Structure Viewer"))
        main_layout.addWidget(self.structure_card())

        self.eq_fig, self.eq_ax = None, None   # Equilibrium concentration plot

    # ------------------------------
    # UI BUILDERS
    # ------------------------------
    def section_label(self, text):
        label = QLabel(text)
        label.setProperty("section", True)
        return label

    def create_card(self, layout_inside):
        card = QWidget()
        card.setLayout(layout_inside)
        #card.setStyleSheet("background-color: white; border-radius: 10px; padding: 10px;")
        card.setObjectName("card")
        return card

    def model_card(self):
        layout = QVBoxLayout()
        layout.setSpacing(2)

        # Material
        row = QHBoxLayout()
        self.dna_radio = QRadioButton("DNA")
        self.rna_radio = QRadioButton("RNA")
        self.dna_radio.setChecked(True)
        row.addWidget(QLabel("Material:"))
        row.addWidget(self.dna_radio)
        row.addWidget(self.rna_radio)
        layout.addLayout(row)

        self.dna_radio.toggled.connect(self.on_material_change)

        # Temperature
        self.temperature = QDoubleSpinBox()
        self.temperature.setValue(37)
        #self.temperature.setSuffix(" °C")
        layout.addLayout(self.form_row("Temperature (ºC):", self.temperature))

        # Sodium
        self.sodium = QDoubleSpinBox()
        self.sodium.setValue(1.0)
        self.sodium.setRange(0.05, 1.10) 
        #self.sodium.setSuffix(" M")
        layout.addLayout(self.form_row("Sum of monovalent Na\u207A / K\u207A / NH\u2084\u207A ions (0.05 - 1.10 M):", self.sodium))


        # Magnesium
        self.magnesium = QDoubleSpinBox()
        self.magnesium.setValue(0.0)
        self.magnesium.setRange(0, 0.2) 
        #self.magnesium.setSuffix(" M")
        layout.addLayout(self.form_row("Sum of divalent Mg\u00b2\u207A ions (0 - 0.2 M):", self.magnesium))

        return self.create_card(layout)

    def strand_card(self):
        layout = QVBoxLayout()

        self.strand_container = QVBoxLayout()
        self.strand_container.setSpacing(2)
        self.strand_rows = []

        self.add_strand()

        #scroll_widget = QWidget()
        #scroll_widget.setLayout(self.strand_container)
        #scroll = QScrollArea()
        #scroll.setWidgetResizable(True)
        #scroll.setWidget(scroll_widget)
        #scroll.setFixedHeight(180)
        #layout.addWidget(scroll)
        
        layout.addLayout(self.strand_container)

        self.add_strand_btn = QPushButton("Add Strand")
        self.add_strand_btn.setProperty("secondary", True)
        self.add_strand_btn.clicked.connect(self.add_strand)

        layout.addWidget(self.add_strand_btn)

        # Max size
        self.max_size = QSpinBox()
        self.max_size.setValue(4)
        #self.max_size.setSuffix(" strands")
        layout.addLayout(self.form_row("Max complex size (strands):", self.max_size))

        # Run button
        self.run_btn = QPushButton("Run Analysis")
        self.run_btn.clicked.connect(self.run_analysis)
        layout.addWidget(self.run_btn)

        return self.create_card(layout)

    def structure_card(self):
        layout = QVBoxLayout()
        layout.setSpacing(2)

        self.complex_dropdown = QComboBox()
        self.complex_dropdown.setFixedWidth(150)
        
        self.base_size = QComboBox()
        self.base_size.addItems(["big", "small"])
        self.base_size.setFixedWidth(80)
        
        self.rotate = QSpinBox()
        self.rotate.setRange(0, 180)
        self.rotate.setValue(80)
        self.rotate.setSuffix("º")

        layout.addLayout(self.form_row("Complex:", self.complex_dropdown))
        layout.addLayout(self.form_row("Base size:", self.base_size))
        layout.addLayout(self.form_row("Rotation:", self.rotate))

        self.plot_btn = QPushButton("Plot Structure")
        self.plot_btn.clicked.connect(self.plot_structure)
        layout.addWidget(self.plot_btn)

        return self.create_card(layout)

    # arrange widgets into a row
    def form_row(self, label_text, widget):
        row = QHBoxLayout()
        row.addWidget(QLabel(label_text))
        row.addStretch()
        row.addWidget(widget)
        return row

    def is_dna(self):
        return self.dna_radio.isChecked()

    # ------------------------------
    # LOGIC
    # ------------------------------
    # add strand button adds new row
    def add_strand(self):
        row = StrandRow(self)
        self.strand_rows.append(row)
        self.strand_container.addWidget(row)

    # if change material, change sequence and model inputs
    def on_material_change(self):
        new_material = "RNA" if self.rna_radio.isChecked() else "DNA"

        for row in self.strand_rows:
            row.convert_sequence(new_material)

        if new_material == "RNA":
            self.sodium.setValue(1.0)
            self.magnesium.setValue(0.0)
            self.sodium.setDisabled(True)
            self.magnesium.setDisabled(True)
        else:
            self.sodium.setDisabled(False)
            self.magnesium.setDisabled(False)

    # build dataframe from user input data
    def build_df(self):
        data = []
        valid_bases = "ATGC" if self.is_dna() else "AUGC"

        for i, row in enumerate(self.strand_rows):
            name = row.name.text() or f"strand{i+1}"
            seq = row.sequence.text()
            conc = row.conc.value()

            valid = all(c in valid_bases for c in seq)
            data.append([name, seq, conc, valid])

        return pd.DataFrame(data, columns=["Name", "Sequence", "Conc", "Valid"])

    # calculate equilibrium concentrations
    def run_analysis(self):
        df = self.build_df()

        if df.empty or not df["Valid"].all():
            plt.figure()
            plt.text(0.5,0.5,"Invalid strand input",ha='center')
            plt.show()
            return

        tube_strands = {Strand(seq, name=name): conc * 1e-6
            for name, seq, conc in zip(df["Name"], df["Sequence"], df["Conc"])}
        # set up tube
        t1 = Tube(strands=tube_strands, complexes=SetSpec(max_size=self.max_size.value()), name="Tube t1")
        # set up model
        model = Model(material = "RNA" if self.rna_radio.isChecked() else "DNA", ensemble="stacking", 
                      celsius=float(self.temperature.value()), sodium=float(self.sodium.value()), magnesium=float(self.magnesium.value()))
        # run analysis
        result = tube_analysis(tubes=[t1], compute=["pairs","mfe"], model=model)

        # filter out minor species
        concs = result[t1].complex_concentrations
        filtered = {k:v for k,v in concs.items() if v>1e-8}
        # clean strand name labels
        def clean_name(c):
            return str(c).replace("<Complex ","").replace(">","").replace("(","").replace(")","")
        # sort species into alphabetical order
        sorted_items = sorted(filtered.items(),
            key=lambda x: (clean_name(x[0]).count("+"), clean_name(x[0])))
        
        labels = [clean_name(k) for k, v in sorted_items]
        values = [float(v)*1e6 for k, v in sorted_items]
        
        # Plot bar chart
        # if first time, create plot. if >1 time, erase and replot in same window
        if self.eq_fig is None:
            self.eq_fig, self.eq_ax = plt.subplots(figsize=(10, max(1.5, len(labels)*0.5)))
        else:
            self.eq_ax.clear()  # clear previous bars
    
        bars = self.eq_ax.barh(labels, values, color='red')
        self.eq_ax.invert_yaxis()
        self.eq_ax.set_xlabel("Equilibrium concentration (µM)")
        self.eq_ax.set_ylabel("Complex")
        self.eq_ax.set_xlim(0, max(values)+0.5)

        # add annotation to end of bars
        for bar, val in zip(bars, values):
            self.eq_ax.text(bar.get_width()+0.01, bar.get_y()+bar.get_height()/2, f"{val:.2g} µM", va="center", ha="left", fontsize=10)
    
        self.eq_fig.tight_layout()
    
        self.eq_fig.canvas.draw_idle()  # redraw the existing figure in place
        self.eq_fig.show()
        
        # STORE CLEAN MAPPING 
        self.tube_result_global = result
        self.complex_map = {}

        for complex_obj, conc in sorted_items:
            name = clean_name(complex_obj)
            self.complex_map[name] = complex_obj
        
        # UPDATE DROPDOWN WITH CLEAN NAMES
        self.complex_dropdown.clear()
        self.complex_dropdown.addItems(list(self.complex_map.keys()))

    # plot structure of chosen complex from equilibrium analysis
    def plot_structure(self):
        if self.tube_result_global is None:
            return

        name = self.complex_dropdown.currentText()
        complex_obj = self.complex_map[name]

        walker = self.tube_result_global[complex_obj]
        mfe = min(walker.mfe, key=lambda x: x.energy)

        # list of all strands from user input
        df = self.build_df()
        strand_dict = {row["Name"]: row["Sequence"] for _,row in df.iterrows()}
        strand_names = name.strip("()").split("+")
        full_seq = "".join([strand_dict[name] for name in strand_names])

        # generate coordinates & plot
        nc = NucDraw(str(mfe.structure))
        nc.generate(degree=self.rotate.value())

        nc.plotter(12, bckwargs={'lw':2,'color':'k'}, bpkwargs={'lw':3,'c':'k'}, scwargs={'s':0.01,'c':'k'})
        title=f'MFE proxy structure for {"RNA" if self.rna_radio.isChecked() else "DNA"} complex {name} at {self.temperature.text()}ºC is {mfe.energy:.2f} kcal/mol:'
        plt.title(title)

        # Define base colour and labels
        base_colors_list = ['green','red','blue','black']
        base_order = ['A','T','C','G']
        if self.rna_radio.isChecked():
            base_order = ['A','U','C','G']
            full_seq = full_seq.replace('T','U')
        base_colors = dict(zip(base_order, base_colors_list))
        
        # change base size
        coords = nc.coords
        if self.base_size.currentText() == 'small':
            base_diam = 30
            font_size = 5
        elif self.base_size.currentText() == 'big':
            base_diam = 80
            font_size = 6

        # draw bases
        for i,(x,y) in enumerate(coords):
            base = full_seq[i]
            plt.scatter(x,y, s=base_diam, facecolor=base_colors[base], edgecolor='none', zorder=3)
            plt.text(x-0.3, y-0.1, base, ha='center', va='center', fontsize=font_size, color='white', zorder=4)
            
        # add base numbering annotations
        nc.numbering_outside(20,4,10,{'fontsize':10,'color':'k'})
        
        # if >1 strand, colour backbones & add legend
        if len(strand_names) > 1:
            colors = ['black','blue','red','green','orange','darkviolet','brown']
            strand_colors = colors[:len(strand_names)]
            nc.multistrand_coloring(clr=strand_colors,bckwargs={'lw':2})
            legend_elements = [Line2D([0],[0],color=strand_colors[i],lw=3,label=strand_names[i])
                for i in range(len(strand_names))]
            plt.legend(handles=legend_elements,title='Strand:',loc='upper left')

        # colour base pairing lines by probability of mfe state
        ax = plt.gca()
        pairs = walker.pairs.to_array()
        dot = str(mfe.structure)
        dot = dot.replace('+','')
        stack = []
        pair_map = {}
        for i,c in enumerate(dot):
            if c == '(':
                stack.append(i)
            elif c == ')':
                j = stack.pop()
                pair_map[i] = j
                pair_map[j] = i
        cmap = cm.jet
        norm = mcolors.Normalize(vmin=0,vmax=1)
        for i,j in pair_map.items():
            if i < j:
                p = pairs[i,j]
                x1,y1 = coords[i]
                x2,y2 = coords[j]
                plt.plot([x1,x2], [y1,y2], color=cmap(norm(p)), linewidth=3, zorder=1)
        sm = cm.ScalarMappable(norm=norm,cmap=cmap)
        plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04, label="Probability of MFE state")
        plt.gcf().set_size_inches(7,4)
        plt.show()
    
    # ------------------------------
    # STYLE
    # ------------------------------
    def get_styles(self):
        return """
        QWidget#card {background-color: white; border-radius: 10px; padding: 10px;}

        QLabel[section="true"] {font-size: 16px; font-weight: bold; color: #2f3640;}

        QLineEdit, QComboBox{background-color: white; border: 1px solid #dcdde1; border-radius: 6px; padding: 6px;}

        QPushButton {background-color: #4078f2; color: white; border-radius: 6px; padding: 8px; font-weight: bold; border: none;}
        
        QPushButton:hover {background-color: #2f5fd0;}
        
        QPushButton:pressed {background-color: #244bb5;}
        
        QPushButton {background-clip: padding;}

        QPushButton[secondary="true"] {background-color: #dcdde1; color: #2f3640;}

        QPushButton[secondary="true"]:hover {background-color: #c8cdd8;}
        
        QPushButton[secondary="true"]:pressed {background-color: #b0b5c0;}

        QPushButton[danger="true"] {background-color: #e84118;}

        QPushButton[danger="true"]:hover {background-color: #c23616;}
        
        QPushButton[danger="true"]:pressed {background-color: #a12a13;}

        QSpinBox, QDoubleSpinBox {background-color: white; border: 1px solid #dcdde1; border-radius: 6px; padding: 4px; min-height: 14px;}

        QSpinBox::up-button, QSpinBox::down-button, QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {width: 0px;}

        QComboBox {border: 1px solid #dcdde1; border-radius: 6px; padding: 4px 6px; color: #2f3640;}
        
        QComboBox QAbstractItemView {selection-background-color: #4078f2; selection-color: white; color: #2f3640;}
        """


# ------------------------------
# RUN
# ------------------------------
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = NupackApp()
    window.show()
    sys.exit(app.exec())