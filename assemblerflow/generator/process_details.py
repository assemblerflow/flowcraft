def colored_print(color_string, msg, end_char="\n"):
    print("\x1b[{}{}\x1b[0m".format(color_string, msg), end=end_char)


def procs_dict_parser(procs_dict):
    colored_print("1;32m", "===== L I S T   O F   P R O C E S S E S =====")
    for template, dict_proc_info in procs_dict.items():
        template_str = "\n=> {}".format(template)
        colored_print("1;36m", template_str)
        for info in dict_proc_info:
            info_str = "\t{}: ".format(info)
            if isinstance(dict_proc_info[info], list):
                if len(dict_proc_info[info]) != 0:
                    colored_print("1;34m", info_str, end_char="")
                    print(", ".join(dict_proc_info[info]))
                    #print("\t{}: {}".format(info, ", ".join(
                    #    dict_proc_info[info])))
            else:
                colored_print("1;34m", info_str, end_char="")
                print(dict_proc_info[info])


def proc_collector(process_map, arguments_list):
    # dict to store only the required entries
    procs_dict = {}
    # loops between all process_map Processes
    for name, cls in process_map.items():
        # instantiates each Process class
        cls_inst = cls(template=name)
        # dictionary comprehension to store only the required keys provided by
        # argument_list
        d = {arg_key: vars(cls_inst)[arg_key] for arg_key in vars(cls_inst)
              if arg_key in arguments_list}
        procs_dict[name] = d

    procs_dict_parser(procs_dict)


