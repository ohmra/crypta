FROM julia:latest
ADD "." "/Mcrypt"
WORKDIR "/Mcrypt"
RUN julia -e 'using Pkg; Pkg.add(["Nemo", "Random", "AlgebraicSolving"]);';
CMD ["julia"]
EXPOSE 3000
