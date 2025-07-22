# Base image with Python and common data science tools
FROM python:3.10-slim

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    git \
    curl \
    wget \
    libglib2.0-0 \
    libsm6 \
    libxrender1 \
    libxext6 \
    ca-certificates \
    libxml2-dev \
    libz-dev \
    libglpk-dev \
    libgmp-dev \
    libigraph-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Install Scanpy and dependencies including igraph and leidenalg
RUN pip install --upgrade pip && \
    pip install scanpy anndata pandas numpy matplotlib seaborn jupyter && \
    pip install python-igraph leidenalg louvain

# Set working directory inside the container
WORKDIR /app