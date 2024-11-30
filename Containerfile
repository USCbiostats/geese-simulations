FROM rocker/geospatial:4.4.0

RUN install2.r aphylo
RUN installGithub.r USCbiostats/geese@master

CMD ["bash"]