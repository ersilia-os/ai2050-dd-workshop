# App Guidelines

This folder contains a day-by-day app-based guided activity. Briefly, each subfolder should contain:
- `app.py`: this file is the actual app to be deployed
- `info.py`: contains information such as model urls and data paths
- `utils.py`: contains specific functions to be called by the app
- `plots.py`: contains plotting functions to ensure standard aspect of the apps across different days
- `/data`: all data files necessaary must be here. Remember to use Git LFS if your files exceed the limit on GitHub repositories.
- `requirements.txt`: please list requirements for your app in this folder.


Please, abide by the following guidelines to create the app:
- Do not write functions in the `app.py` file. Only general functions, like the one calling models, should be instantiated at the top. Any other needed functions must go into `utils.py`or `plots.py`.
- Streamlit caching: using the caching properties of streamlit will significantly improve app performance. Make sure to use the decorators appropriately and double test that they are working (i.e, data is cached when expected, and refreshed when needed).
- Session state: the function session_state is used to store information for the current run. Please use it to record button status, model results etc in combination with the caching. An example is provided.
- Streamlit columns: use the columns functionality to create the appropriate spacing. The app properties are set to wide, a maximum of 4 columns is allowed.
- When creating new plots, be sure to use the same style and color coding as the already developed plots.
- Use click buttons to progress to the next section, so that the app appears interactively to the participants. An example is provided.
- Keep the data under the data folder and call it using the dictionary specified in `info.py`
- If any of the files or functions existing in the templates (for example the `example.csv` file) is not being used, please delete them to avoid confusion.
- Do not modify the Dockerfile placed in the root of the folder. Add as needed packages in the `requirements.txt` file
- Please keep in mind that no file outside the subfolder will be read during app deployment.
- The .streamlit folder in the root of the repository is unique to all deployed apps. Do not modify it.
- The .do folder in the root of the repository is unique to all deployed apps. Do not modify it.
