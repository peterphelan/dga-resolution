FROM sagemath/sagemath:latest

WORKDIR /workspace

COPY *.py ./

USER root
RUN chown -R sage:sage /workspace
USER sage

CMD ["sage", "main.py"]
