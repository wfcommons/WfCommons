#!/usr/bin/env bash

_qiime_completion()
{
  local COMP_WORDS=(${COMP_WORDS[*]})
  local incomplete
  if [[ ${COMP_CWORD} -lt 0 ]] ; then
    COMP_CWORD="${#COMP_WORDS[*]}"
    incomplete=""
  else
    incomplete="${COMP_WORDS[COMP_CWORD]}"
  fi

  local curpos nextpos nextword
  nextpos=0

  curpos=${nextpos}
  while :
  do
    nextpos=$((curpos + 1))
    nextword="${COMP_WORDS[nextpos]}"
    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
      if [[ ${incomplete} == -* ]] ; then
        echo "$(compgen -W "--help --version" -- $incomplete)"
      else
        echo "$(compgen -W "info tools dev alignment composition cutadapt dada2 deblur demux diversity emperor feature-classifier feature-table fragment-insertion gneiss longitudinal metadata phylogeny quality-control quality-filter sample-classifier taxa vsearch" -- $incomplete)"
      fi
      return 0
    else
      case "${nextword}" in
        info)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help" -- $incomplete)"
              else
                echo "$(compgen -W "" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        tools)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help" -- $incomplete)"
              else
                echo "$(compgen -W "citations export extract import inspect-metadata peek validate view" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                citations)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                export)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --input-path --output-path --output-format" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                extract)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --input-path --output-path" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                import)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --type --input-path --output-path --input-format --show-importable-types --show-importable-formats" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                inspect-metadata)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --tsv --no-tsv" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                peek)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                validate)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --level" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                view)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --index-extension" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        dev)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help" -- $incomplete)"
              else
                echo "$(compgen -W "export-default-theme import-theme refresh-cache reset-theme" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                export-default-theme)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --output-path" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                import-theme)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --theme" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                refresh-cache)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                reset-theme)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --yes" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        alignment)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "mafft mask" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                mafft)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --p-n-threads --p-parttree --p-no-parttree --o-alignment --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                mask)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-alignment --p-max-gap-frequency --p-min-conservation --o-masked-alignment --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        composition)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "add-pseudocount ancom" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                add-pseudocount)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-pseudocount --o-composition-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                ancom)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --m-metadata-column --p-transform-function --p-difference-function --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        cutadapt)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "demux-paired demux-single trim-paired trim-single" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                demux-paired)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-seqs --m-forward-barcodes-file --m-forward-barcodes-column --m-reverse-barcodes-file --m-reverse-barcodes-column --p-error-rate --p-batch-size --p-minimum-length --o-per-sample-sequences --o-untrimmed-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                demux-single)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-seqs --m-barcodes-file --m-barcodes-column --p-error-rate --p-batch-size --p-minimum-length --o-per-sample-sequences --o-untrimmed-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                trim-paired)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demultiplexed-sequences --p-cores --p-adapter-f --p-front-f --p-anywhere-f --p-adapter-r --p-front-r --p-anywhere-r --p-error-rate --p-indels --p-no-indels --p-times --p-overlap --p-match-read-wildcards --p-no-match-read-wildcards --p-match-adapter-wildcards --p-no-match-adapter-wildcards --p-minimum-length --p-discard-untrimmed --p-no-discard-untrimmed --o-trimmed-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                trim-single)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demultiplexed-sequences --p-cores --p-adapter --p-front --p-anywhere --p-error-rate --p-indels --p-no-indels --p-times --p-overlap --p-match-read-wildcards --p-no-match-read-wildcards --p-match-adapter-wildcards --p-no-match-adapter-wildcards --p-minimum-length --p-discard-untrimmed --p-no-discard-untrimmed --o-trimmed-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        dada2)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "denoise-paired denoise-pyro denoise-single" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                denoise-paired)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demultiplexed-seqs --p-trunc-len-f --p-trunc-len-r --p-trim-left-f --p-trim-left-r --p-max-ee-f --p-max-ee-r --p-trunc-q --p-chimera-method --p-min-fold-parent-over-abundance --p-n-threads --p-n-reads-learn --p-hashed-feature-ids --p-no-hashed-feature-ids --o-table --o-representative-sequences --o-denoising-stats --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                denoise-pyro)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demultiplexed-seqs --p-trunc-len --p-trim-left --p-max-ee --p-trunc-q --p-max-len --p-chimera-method --p-min-fold-parent-over-abundance --p-n-threads --p-n-reads-learn --p-hashed-feature-ids --p-no-hashed-feature-ids --o-table --o-representative-sequences --o-denoising-stats --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                denoise-single)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demultiplexed-seqs --p-trunc-len --p-trim-left --p-max-ee --p-trunc-q --p-chimera-method --p-min-fold-parent-over-abundance --p-n-threads --p-n-reads-learn --p-hashed-feature-ids --p-no-hashed-feature-ids --o-table --o-representative-sequences --o-denoising-stats --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        deblur)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "denoise-16S denoise-other visualize-stats" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                denoise-16S)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demultiplexed-seqs --p-trim-length --p-left-trim-len --p-sample-stats --p-no-sample-stats --p-mean-error --p-indel-prob --p-indel-max --p-min-reads --p-min-size --p-jobs-to-start --p-hashed-feature-ids --p-no-hashed-feature-ids --o-table --o-representative-sequences --o-stats --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                denoise-other)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demultiplexed-seqs --i-reference-seqs --p-trim-length --p-left-trim-len --p-sample-stats --p-no-sample-stats --p-mean-error --p-indel-prob --p-indel-max --p-min-reads --p-min-size --p-jobs-to-start --p-hashed-feature-ids --p-no-hashed-feature-ids --o-table --o-representative-sequences --o-stats --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                visualize-stats)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-deblur-stats --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        demux)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "emp-paired emp-single filter-samples subsample-paired subsample-single summarize" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                emp-paired)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-seqs --m-barcodes-file --m-barcodes-column --p-golay-error-correction --p-no-golay-error-correction --p-rev-comp-barcodes --p-no-rev-comp-barcodes --p-rev-comp-mapping-barcodes --p-no-rev-comp-mapping-barcodes --o-per-sample-sequences --o-error-correction-details --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                emp-single)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-seqs --m-barcodes-file --m-barcodes-column --p-golay-error-correction --p-no-golay-error-correction --p-rev-comp-barcodes --p-no-rev-comp-barcodes --p-rev-comp-mapping-barcodes --p-no-rev-comp-mapping-barcodes --o-per-sample-sequences --o-error-correction-details --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                filter-samples)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demux --m-metadata-file --p-where --p-exclude-ids --p-no-exclude-ids --o-filtered-demux --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                subsample-paired)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --p-fraction --o-subsampled-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                subsample-single)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --p-fraction --o-subsampled-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                summarize)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-data --p-n --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        diversity)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "adonis alpha alpha-correlation alpha-group-significance alpha-phylogenetic alpha-phylogenetic-alt alpha-rarefaction beta beta-correlation beta-group-significance beta-phylogenetic beta-rarefaction bioenv core-metrics core-metrics-phylogenetic filter-distance-matrix mantel pcoa pcoa-biplot procrustes-analysis" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                adonis)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-distance-matrix --m-metadata-file --p-formula --p-permutations --p-n-jobs --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                alpha)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-metric --o-alpha-diversity --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                alpha-correlation)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-alpha-diversity --m-metadata-file --p-method --p-intersect-ids --p-no-intersect-ids --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                alpha-group-significance)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-alpha-diversity --m-metadata-file --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                alpha-phylogenetic)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-phylogeny --p-metric --o-alpha-diversity --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                alpha-phylogenetic-alt)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-phylogeny --p-metric --o-alpha-diversity --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                alpha-rarefaction)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-phylogeny --p-max-depth --p-metrics --m-metadata-file --p-min-depth --p-steps --p-iterations --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                beta)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-metric --p-pseudocount --p-n-jobs --o-distance-matrix --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                beta-correlation)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-distance-matrix --m-metadata-file --m-metadata-column --p-method --p-permutations --p-intersect-ids --p-no-intersect-ids --p-label1 --p-label2 --o-metadata-distance-matrix --o-mantel-scatter-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                beta-group-significance)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-distance-matrix --m-metadata-file --m-metadata-column --p-method --p-pairwise --p-no-pairwise --p-permutations --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                beta-phylogenetic)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-phylogeny --p-metric --p-n-jobs --p-variance-adjusted --p-no-variance-adjusted --p-alpha --p-bypass-tips --p-no-bypass-tips --o-distance-matrix --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                beta-rarefaction)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-phylogeny --p-metric --p-clustering-method --m-metadata-file --p-sampling-depth --p-iterations --p-correlation-method --p-color-scheme --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                bioenv)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-distance-matrix --m-metadata-file --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                core-metrics)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-sampling-depth --m-metadata-file --p-with-replacement --p-no-with-replacement --p-n-jobs --o-rarefied-table --o-observed-otus-vector --o-shannon-vector --o-evenness-vector --o-jaccard-distance-matrix --o-bray-curtis-distance-matrix --o-jaccard-pcoa-results --o-bray-curtis-pcoa-results --o-jaccard-emperor --o-bray-curtis-emperor --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                core-metrics-phylogenetic)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-phylogeny --p-sampling-depth --m-metadata-file --p-n-jobs --o-rarefied-table --o-faith-pd-vector --o-observed-otus-vector --o-shannon-vector --o-evenness-vector --o-unweighted-unifrac-distance-matrix --o-weighted-unifrac-distance-matrix --o-jaccard-distance-matrix --o-bray-curtis-distance-matrix --o-unweighted-unifrac-pcoa-results --o-weighted-unifrac-pcoa-results --o-jaccard-pcoa-results --o-bray-curtis-pcoa-results --o-unweighted-unifrac-emperor --o-weighted-unifrac-emperor --o-jaccard-emperor --o-bray-curtis-emperor --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                filter-distance-matrix)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-distance-matrix --m-metadata-file --p-where --p-exclude-ids --p-no-exclude-ids --o-filtered-distance-matrix --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                mantel)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-dm1 --i-dm2 --p-method --p-permutations --p-intersect-ids --p-no-intersect-ids --p-label1 --p-label2 --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                pcoa)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-distance-matrix --p-number-of-dimensions --o-pcoa --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                pcoa-biplot)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-pcoa --i-features --o-biplot --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                procrustes-analysis)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-reference --i-other --p-dimensions --o-transformed-reference --o-transformed-other --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        emperor)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "biplot plot procrustes-plot" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                biplot)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-biplot --m-sample-metadata-file --m-feature-metadata-file --p-ignore-missing-samples --p-no-ignore-missing-samples --p-number-of-features --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                plot)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-pcoa --m-metadata-file --p-custom-axes --p-ignore-missing-samples --p-no-ignore-missing-samples --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                procrustes-plot)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-reference-pcoa --i-other-pcoa --m-metadata-file --p-custom-axes --p-ignore-missing-samples --p-no-ignore-missing-samples --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        feature-classifier)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "classify-consensus-blast classify-consensus-vsearch classify-hybrid-vsearch-sklearn classify-sklearn extract-reads fit-classifier-naive-bayes fit-classifier-sklearn" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                classify-consensus-blast)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-query --i-reference-reads --i-reference-taxonomy --p-maxaccepts --p-perc-identity --p-query-cov --p-strand --p-evalue --p-min-consensus --p-unassignable-label --o-classification --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                classify-consensus-vsearch)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-query --i-reference-reads --i-reference-taxonomy --p-maxaccepts --p-perc-identity --p-query-cov --p-strand --p-min-consensus --p-unassignable-label --p-search-exact --p-no-search-exact --p-top-hits-only --p-no-top-hits-only --p-threads --o-classification --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                classify-hybrid-vsearch-sklearn)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-query --i-reference-reads --i-reference-taxonomy --i-classifier --p-maxaccepts --p-perc-identity --p-query-cov --p-strand --p-min-consensus --p-reads-per-batch --p-confidence --p-read-orientation --p-threads --p-prefilter --p-no-prefilter --p-sample-size --p-randseed --o-classification --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                classify-sklearn)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-reads --i-classifier --p-reads-per-batch --p-n-jobs --p-pre-dispatch --p-confidence --p-read-orientation --o-classification --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                extract-reads)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --p-f-primer --p-r-primer --p-trunc-len --p-trim-left --p-identity --p-min-length --p-max-length --p-n-jobs --p-batch-size --o-reads --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                fit-classifier-naive-bayes)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-reference-reads --i-reference-taxonomy --i-class-weight --p-classify--alpha --p-classify--chunk-size --p-classify--class-prior --p-classify--fit-prior --p-no-classify--fit-prior --p-feat-ext--alternate-sign --p-no-feat-ext--alternate-sign --p-feat-ext--analyzer --p-feat-ext--binary --p-no-feat-ext--binary --p-feat-ext--decode-error --p-feat-ext--encoding --p-feat-ext--input --p-feat-ext--lowercase --p-no-feat-ext--lowercase --p-feat-ext--n-features --p-feat-ext--ngram-range --p-feat-ext--norm --p-feat-ext--preprocessor --p-feat-ext--stop-words --p-feat-ext--strip-accents --p-feat-ext--token-pattern --p-feat-ext--tokenizer --p-verbose --p-no-verbose --o-classifier --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                fit-classifier-sklearn)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-reference-reads --i-reference-taxonomy --i-class-weight --p-classifier-specification --o-classifier --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        feature-table)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "core-features filter-features filter-samples filter-seqs group heatmap merge merge-seqs merge-taxa presence-absence rarefy relative-frequency subsample summarize tabulate-seqs transpose" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                core-features)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-min-fraction --p-max-fraction --p-steps --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                filter-features)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-min-frequency --p-max-frequency --p-min-samples --p-max-samples --m-metadata-file --p-where --p-exclude-ids --p-no-exclude-ids --o-filtered-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                filter-samples)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-min-frequency --p-max-frequency --p-min-features --p-max-features --m-metadata-file --p-where --p-exclude-ids --p-no-exclude-ids --o-filtered-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                filter-seqs)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-data --i-table --m-metadata-file --p-where --p-exclude-ids --p-no-exclude-ids --o-filtered-data --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                group)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-axis --m-metadata-file --m-metadata-column --p-mode --o-grouped-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                heatmap)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-sample-metadata-file --m-sample-metadata-column --m-feature-metadata-file --m-feature-metadata-column --p-normalize --p-no-normalize --p-title --p-metric --p-method --p-cluster --p-color-scheme --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                merge)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-tables --p-overlap-method --o-merged-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                merge-seqs)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-data --o-merged-data --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                merge-taxa)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-data --o-merged-data --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                presence-absence)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --o-presence-absence-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                rarefy)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-sampling-depth --p-with-replacement --p-no-with-replacement --o-rarefied-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                relative-frequency)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --o-relative-frequency-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                subsample)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-subsampling-depth --p-axis --o-sampled-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                summarize)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-sample-metadata-file --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                tabulate-seqs)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-data --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                transpose)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --o-transposed-feature-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        fragment-insertion)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "classify-otus-experimental filter-features sepp" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                classify-otus-experimental)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-representative-sequences --i-tree --i-reference-taxonomy --o-classification --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                filter-features)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-tree --o-filtered-table --o-removed-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                sepp)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-representative-sequences --i-reference-database --p-alignment-subset-size --p-placement-subset-size --p-threads --p-debug --p-no-debug --o-tree --o-placements --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        gneiss)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "assign-ids balance-taxonomy correlation-clustering dendrogram-heatmap gradient-clustering ilr-hierarchical ilr-phylogenetic lme-regression ols-regression" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                assign-ids)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-input-table --i-input-tree --o-output-table --o-output-tree --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                balance-taxonomy)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-tree --i-taxonomy --p-balance-name --p-pseudocount --p-taxa-level --p-n-features --p-threshold --m-metadata-file --m-metadata-column --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                correlation-clustering)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --p-pseudocount --o-clustering --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                dendrogram-heatmap)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-tree --m-metadata-file --m-metadata-column --p-pseudocount --p-ndim --p-method --p-color-map --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                gradient-clustering)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-gradient-file --m-gradient-column --p-ignore-missing-samples --p-no-ignore-missing-samples --p-weighted --p-no-weighted --o-clustering --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                ilr-hierarchical)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-tree --p-pseudocount --o-balances --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                ilr-phylogenetic)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-tree --p-pseudocount --o-balances --o-hierarchy --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                lme-regression)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-tree --m-metadata-file --p-formula --p-groups --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                ols-regression)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-tree --m-metadata-file --p-formula --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        longitudinal)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "anova feature-volatility first-differences first-distances linear-mixed-effects maturity-index nmit pairwise-differences pairwise-distances plot-feature-volatility volatility" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                anova)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --m-metadata-file --p-formula --p-sstype --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                feature-volatility)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --p-state-column --p-individual-id-column --p-cv --p-random-state --p-n-jobs --p-n-estimators --p-estimator --p-parameter-tuning --p-no-parameter-tuning --p-missing-samples --p-importance-threshold --p-feature-count --o-filtered-table --o-feature-importance --o-volatility-plot --o-accuracy-results --o-sample-estimator --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                first-differences)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --p-state-column --p-individual-id-column --p-metric --p-replicate-handling --p-baseline --o-first-differences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                first-distances)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-distance-matrix --m-metadata-file --p-state-column --p-individual-id-column --p-baseline --p-replicate-handling --o-first-distances --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                linear-mixed-effects)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --p-state-column --p-individual-id-column --p-metric --p-group-columns --p-random-effects --p-palette --p-lowess --p-no-lowess --p-ci --p-formula --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                maturity-index)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --p-state-column --p-group-by --p-control --p-individual-id-column --p-estimator --p-n-estimators --p-test-size --p-step --p-cv --p-random-state --p-n-jobs --p-parameter-tuning --p-no-parameter-tuning --p-optimize-feature-selection --p-no-optimize-feature-selection --p-stratify --p-no-stratify --p-missing-samples --p-feature-count --o-sample-estimator --o-feature-importance --o-predictions --o-model-summary --o-accuracy-results --o-maz-scores --o-clustermap --o-volatility-plots --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                nmit)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --p-individual-id-column --p-corr-method --p-dist-method --o-distance-matrix --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                pairwise-differences)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --p-metric --p-state-column --p-state-1 --p-state-2 --p-individual-id-column --p-group-column --p-parametric --p-no-parametric --p-palette --p-replicate-handling --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                pairwise-distances)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-distance-matrix --m-metadata-file --p-group-column --p-state-column --p-state-1 --p-state-2 --p-individual-id-column --p-parametric --p-no-parametric --p-palette --p-replicate-handling --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                plot-feature-volatility)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-importances --m-metadata-file --p-state-column --p-individual-id-column --p-default-group-column --p-yscale --p-importance-threshold --p-feature-count --p-missing-samples --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                volatility)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --p-state-column --p-individual-id-column --p-default-group-column --p-default-metric --p-yscale --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        metadata)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "distance-matrix tabulate" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                distance-matrix)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --m-metadata-file --m-metadata-column --o-distance-matrix --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                tabulate)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --m-input-file --p-page-size --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        phylogeny)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "align-to-tree-mafft-fasttree fasttree filter-table iqtree iqtree-ultrafast-bootstrap midpoint-root raxml raxml-rapid-bootstrap" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                align-to-tree-mafft-fasttree)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --p-n-threads --p-mask-max-gap-frequency --p-mask-min-conservation --o-alignment --o-masked-alignment --o-tree --o-rooted-tree --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                fasttree)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-alignment --p-n-threads --o-tree --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                filter-table)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-tree --o-filtered-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                iqtree)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-alignment --p-seed --p-n-cores --p-n-runs --p-substitution-model --p-n-init-pars-trees --p-n-top-init-trees --p-n-best-retain-trees --p-n-iter --p-stop-iter --p-perturb-nni-strength --p-spr-radius --p-allnni --p-no-allnni --p-fast --p-no-fast --p-alrt --p-abayes --p-no-abayes --p-lbp --p-safe --p-no-safe --o-tree --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                iqtree-ultrafast-bootstrap)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-alignment --p-seed --p-n-cores --p-n-runs --p-substitution-model --p-bootstrap-replicates --p-n-init-pars-trees --p-n-top-init-trees --p-n-best-retain-trees --p-stop-iter --p-perturb-nni-strength --p-spr-radius --p-n-max-ufboot-iter --p-n-ufboot-steps --p-min-cor-ufboot --p-ep-break-ufboot --p-allnni --p-no-allnni --p-alrt --p-abayes --p-no-abayes --p-lbp --p-bnni --p-no-bnni --p-safe --p-no-safe --o-tree --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                midpoint-root)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-tree --o-rooted-tree --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                raxml)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-alignment --p-seed --p-n-searches --p-n-threads --p-raxml-version --p-substitution-model --o-tree --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                raxml-rapid-bootstrap)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-alignment --p-seed --p-rapid-bootstrap-seed --p-bootstrap-replicates --p-n-threads --p-raxml-version --p-substitution-model --o-tree --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        quality-control)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "evaluate-composition evaluate-seqs evaluate-taxonomy exclude-seqs" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                evaluate-composition)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-expected-features --i-observed-features --p-depth --p-palette --p-plot-tar --p-no-plot-tar --p-plot-tdr --p-no-plot-tdr --p-plot-r-value --p-no-plot-r-value --p-plot-r-squared --p-no-plot-r-squared --p-plot-bray-curtis --p-no-plot-bray-curtis --p-plot-jaccard --p-no-plot-jaccard --p-plot-observed-features --p-no-plot-observed-features --p-plot-observed-features-ratio --p-no-plot-observed-features-ratio --m-metadata-file --m-metadata-column --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                evaluate-seqs)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-query-sequences --i-reference-sequences --p-show-alignments --p-no-show-alignments --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                evaluate-taxonomy)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-expected-taxa --i-observed-taxa --i-feature-table --p-depth --p-palette --p-require-exp-ids --p-no-require-exp-ids --p-require-obs-ids --p-no-require-obs-ids --p-sample-id --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                exclude-seqs)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-query-sequences --i-reference-sequences --p-method --p-perc-identity --p-evalue --p-perc-query-aligned --p-threads --o-sequence-hits --o-sequence-misses --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        quality-filter)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "q-score q-score-joined" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                q-score)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demux --p-min-quality --p-quality-window --p-min-length-fraction --p-max-ambiguous --o-filtered-sequences --o-filter-stats --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                q-score-joined)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demux --p-min-quality --p-quality-window --p-min-length-fraction --p-max-ambiguous --o-filtered-sequences --o-filter-stats --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        sample-classifier)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "classify-samples classify-samples-from-dist classify-samples-ncv confusion-matrix fit-classifier fit-regressor heatmap metatable predict-classification predict-regression regress-samples regress-samples-ncv scatterplot split-table summarize" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                classify-samples)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --m-metadata-column --p-test-size --p-step --p-cv --p-random-state --p-n-jobs --p-n-estimators --p-estimator --p-optimize-feature-selection --p-no-optimize-feature-selection --p-parameter-tuning --p-no-parameter-tuning --p-palette --p-missing-samples --o-sample-estimator --o-feature-importance --o-predictions --o-model-summary --o-accuracy-results --o-probabilities --o-heatmap --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                classify-samples-from-dist)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-distance-matrix --m-metadata-file --m-metadata-column --p-k --p-palette --o-predictions --o-accuracy-results --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                classify-samples-ncv)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --m-metadata-column --p-cv --p-random-state --p-n-jobs --p-n-estimators --p-estimator --p-parameter-tuning --p-no-parameter-tuning --p-missing-samples --o-predictions --o-feature-importance --o-probabilities --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                confusion-matrix)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-predictions --i-probabilities --m-truth-file --m-truth-column --p-missing-samples --p-vmin --p-vmax --p-palette --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                fit-classifier)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --m-metadata-column --p-step --p-cv --p-random-state --p-n-jobs --p-n-estimators --p-estimator --p-optimize-feature-selection --p-no-optimize-feature-selection --p-parameter-tuning --p-no-parameter-tuning --p-missing-samples --o-sample-estimator --o-feature-importance --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                fit-regressor)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --m-metadata-column --p-step --p-cv --p-random-state --p-n-jobs --p-n-estimators --p-estimator --p-optimize-feature-selection --p-no-optimize-feature-selection --p-parameter-tuning --p-no-parameter-tuning --p-missing-samples --o-sample-estimator --o-feature-importance --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                heatmap)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-importance --m-sample-metadata-file --m-sample-metadata-column --m-feature-metadata-file --m-feature-metadata-column --p-feature-count --p-importance-threshold --p-group-samples --p-no-group-samples --p-normalize --p-no-normalize --p-missing-samples --p-metric --p-method --p-cluster --p-color-scheme --o-heatmap --o-filtered-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                metatable)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --p-missing-samples --p-missing-values --p-drop-all-unique --p-no-drop-all-unique --o-converted-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                predict-classification)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-sample-estimator --p-n-jobs --o-predictions --o-probabilities --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                predict-regression)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-sample-estimator --p-n-jobs --o-predictions --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                regress-samples)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --m-metadata-column --p-test-size --p-step --p-cv --p-random-state --p-n-jobs --p-n-estimators --p-estimator --p-optimize-feature-selection --p-no-optimize-feature-selection --p-stratify --p-no-stratify --p-parameter-tuning --p-no-parameter-tuning --p-missing-samples --o-sample-estimator --o-feature-importance --o-predictions --o-model-summary --o-accuracy-results --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                regress-samples-ncv)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --m-metadata-column --p-cv --p-random-state --p-n-jobs --p-n-estimators --p-estimator --p-stratify --p-no-stratify --p-parameter-tuning --p-no-parameter-tuning --p-missing-samples --o-predictions --o-feature-importance --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                scatterplot)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-predictions --m-truth-file --m-truth-column --p-missing-samples --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                split-table)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --m-metadata-file --m-metadata-column --p-test-size --p-random-state --p-stratify --p-no-stratify --p-missing-samples --o-training-table --o-test-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                summarize)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sample-estimator --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        taxa)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "barplot collapse filter-seqs filter-table" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                barplot)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-taxonomy --m-metadata-file --o-visualization --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                collapse)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-taxonomy --p-level --o-collapsed-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                filter-seqs)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --i-taxonomy --p-include --p-exclude --p-query-delimiter --p-mode --o-filtered-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                filter-table)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-table --i-taxonomy --p-include --p-exclude --p-query-delimiter --p-mode --o-filtered-table --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

        vsearch)
          curpos=${nextpos}
          while :
          do
            nextpos=$((curpos + 1))
            nextword="${COMP_WORDS[nextpos]}"
            if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
              if [[ ${incomplete} == -* ]] ; then
                echo "$(compgen -W "--help --version --citations" -- $incomplete)"
              else
                echo "$(compgen -W "cluster-features-closed-reference cluster-features-de-novo cluster-features-open-reference dereplicate-sequences join-pairs uchime-denovo uchime-ref" -- $incomplete)"
              fi
              return 0
            else
              case "${nextword}" in
                cluster-features-closed-reference)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --i-table --i-reference-sequences --p-perc-identity --p-strand --p-threads --o-clustered-table --o-clustered-sequences --o-unmatched-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                cluster-features-de-novo)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --i-table --p-perc-identity --p-threads --o-clustered-table --o-clustered-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                cluster-features-open-reference)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --i-table --i-reference-sequences --p-perc-identity --p-strand --p-threads --o-clustered-table --o-clustered-sequences --o-new-reference-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                dereplicate-sequences)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --p-derep-prefix --p-no-derep-prefix --o-dereplicated-table --o-dereplicated-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                join-pairs)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-demultiplexed-seqs --p-truncqual --p-minlen --p-maxns --p-allowmergestagger --p-no-allowmergestagger --p-minovlen --p-maxdiffs --p-minmergelen --p-maxmergelen --p-maxee --p-qmin --p-qminout --p-qmax --p-qmaxout --p-threads --o-joined-sequences --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                uchime-denovo)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --i-table --p-dn --p-mindiffs --p-mindiv --p-minh --p-xn --o-chimeras --o-nonchimeras --o-stats --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

                uchime-ref)
                  curpos=${nextpos}
                  while :
                  do
                    nextpos=$((curpos + 1))
                    nextword="${COMP_WORDS[nextpos]}"
                    if [[ ${nextpos} -eq ${COMP_CWORD} ]] ; then
                      if [[ ${incomplete} == -* ]] ; then
                        echo "$(compgen -W "--help --i-sequences --i-table --i-reference-sequences --p-dn --p-mindiffs --p-mindiv --p-minh --p-xn --p-threads --o-chimeras --o-nonchimeras --o-stats --output-dir --verbose --quiet --citations" -- $incomplete)"
                      else
                        echo "$(compgen -W "" -- $incomplete)"
                      fi
                      return 0
                    else
                      case "${nextword}" in

                      esac
                      curpos=${nextpos}
                    fi
                  done
                  ;;

              esac
              curpos=${nextpos}
            fi
          done
          ;;

      esac
      curpos=${nextpos}
    fi
  done

  return 0
}

_qiime_completion
