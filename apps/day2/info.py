about = [
    "This app is part of the Schmidt Futures AI2050 program.",
    "The workshop been jointly developed by the [Ersilia Open Source Initiative](https://ersilia.io) and the [H3D Foundation](https://h3dfoundation.org/).",
    "For more information about this workshop, please see code and data in this [GitHub repository](https://github.com/ersilia-os/ai2050-dd-workshop).",
    "If you have a more advanced dataset in mind or a use case for your research, please contact us at: [hello@ersilia.io](mailto:hello@ersilia.io)."
]

intro = """
Today we will focus on understanding and practicing the basic steps to train an ML model. We will be given the results of an experiment for _A.baumannii_ inhibition (growth), and based on these we will test different options (classifiers, regressors) as well as assessing different molecular descriptors.
""".strip()

exp  = """You are handed, directly from your experimental collaborators, the following dataset:"""

model_urls = {
    "eos4tw0": ["https://eos9ei3-tkreo.ondigitalocean.app/"
                ],
                
    "eos4u6p": ["https://eos43at-zqx9x.ondigitalocean.app/"
                ],
    "eos3804": ["https://eos3894-gz5nz.ondigitalocean.app/"]
}

library_filenames = {
    "Compound library 1": "abaumannii_subset250_0.csv",
    "Compound library 2": "abaumannii_subset250_1.csv",
    "Compound library 3": "abaumannii_subset250_2.csv",
    "Compound library 4": "abaumannii_subset250_3.csv",
}

q1 = [
    "- What does each row represent?",
    "- Which are the columns with experimental values?",
    "- There is a vital piece of information missing in this table, what is it?",
    "- What is the experimental assay measuring?",
    "- Do we want to obtain higher or lower experimental values?",
    "- Which columns are we going to use moving forward?"
]

q2 = [
    "- Why do we need to decide an activity cut-off?",
    "- What does 0 indicate? And 1?",
    "- Why to we want to know the Mean OD of the dataset?",
    "- What else do we need to take into account to define a cut-off?",
    "- What would be a good cut-off in this case?",
    "- Is it a balanced dataset? Why or why not?",
    "- What is the author's defined activity cut-off?",
]

q3 = [
    "- What do 1D, 2D and 3D descriptors take into account?",
    "- What do Morgan Fingerprints represent?",
    "- What do Chemical Checker signatures represent?",
    "- Which one would be most appropriate for this case?",
]

q4 = [
    "- What is a cross-validation experiment?"
    "- Why do we need a train/test split?",
    "- What could be a good ratio for the train/test split?",
    "- What do we have to watch out for in this particular dataset?",
    "- Which metrics can we use to evaluate model performance?",
]

q4_followup = [
    "- Are these good performances?",
    "- Which descriptors seem to work better?",
    "- What would happen if we increase or decrease the cut-off?",
    "- How does our quick modelling compare to the author's work?",
    "- What else could we do now?"
]

q5 = [
    "- What is the output of a classifier?",
    "- What do we consider an Active in the predictions?",
    "- What did the authors consider an Active in their predictions?",
    "- Are these results expected?"
]