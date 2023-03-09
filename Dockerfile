FROM amrit3701/freecad-cli:0.20-amd64

COPY requirements.txt .
RUN python3.8 -m pip install -r requirements.txt

WORKDIR /app
COPY . /app

# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
USER appuser

# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
CMD ["python3.8", "puzzlePiecesGenerator.py"]
