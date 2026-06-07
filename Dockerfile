FROM python:3.12-slim

WORKDIR /app
COPY . .
RUN pip install --no-cache-dir .

# Users bind-mount their data under /data
WORKDIR /data
ENTRYPOINT ["scagg"]
