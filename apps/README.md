# App Guidelines

This folder contains a day-by-day app-based guided activity. Briefly, each folder should contain:
- `app.py`: this file is the actual app to be deployed
- `info.p`: contains information such as model urls and data paths
- `utils.py`: contains specific functions to be called by the app
- `plots.py`: contains plotting functions to ensure standard aspect of the apps across different days


Please, abide by the following guidelines to create the app:
- Do not write functions in the `app.py` file. Only general functions, like the one calling models, should be instantiated at the top. Any other needed functions must go into `utils.py`or `plots.py`.
- When creating new plots, be sure to use the same style and color coding as the already developed plots.
- Streamlit caching: using the caching properties of streamlit will significantly improve app performance. Make sure to use the decorators appropriately and double test that they are working (i.e, data is cached when expected, and refreshed when needed).
- Streamlit columns: use the columns functionality to create the appropriate spacing. The app properties are set to wide, a maximum of 4 columns is allowed.
- Use click buttons to progress to the next section, so that the app appears interactively to the participants. An example is provided
- Keep the data under data/dayX folder and call it using a dictionary specified in info.py
- If any of the files or functions existing in the template folder is not being used, please delete them to avoid confusion.