FROM python:3.10.7-bullseye

WORKDIR .

COPY . .

RUN python -m pip install --upgrade pip
RUN python -m pip install -r requirements.txt
RUN apt-get update && \
    apt-get install -y git git-lfs && \
    git lfs install
RUN git lfs pull

CMD ["streamlit", "run", "app.py"]