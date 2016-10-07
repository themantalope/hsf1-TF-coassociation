import pandas as pd
import os

data_dir = os.path.join("..", "data_files")
expression_file = os.path.join(data_dir, "hsf1_reg_oscillatory_genes.xlsx")
cc_tf_binding_file = os.path.join(data_dir, "fischer_decaprio_cell_cycle_TFs_targets.xlsx")

osc_expression = pd.read_excel(open(expression_file))
cc_tf_binding = pd.read_excel(open(cc_tf_binding_file), sheetname=None)


for sheet in cc_tf_binding.keys():
    cursheet = cc_tf_binding[sheet]
    tf_name = sheet.replace("targets", "")
    tf_name = tf_name.strip()

    target_ensembl_ids = cursheet["ensembl ID"].values.tolist()

    osc_expression[tf_name] = osc_expression.apply(lambda x: 1 if x["ENSG"] in target_ensembl_ids else 0, axis=1)

writer = pd.ExcelWriter(os.path.join(data_dir, "hsf1_reg_oscillatory_genes.cctf.xlsx"))
osc_expression.to_excel(writer)
writer.save()
