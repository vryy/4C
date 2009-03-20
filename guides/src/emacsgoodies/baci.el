;; baci.el
;; $Id: baci.el,v 1.1 2006/11/20 01:16:04 pbb Exp $
;;
;; About: A simple Emacs mode for Baci input DAT files
;;
;; Typist: B Bornemann
;; Notice: File comes without any warranty what so ever.
;;
;; Based on: Jorgen Bergstrom's lsdyna.el <jorgen@polymerFEM.com>
;; URL:      http://polymerFEM.com
;; 
;; Installation:
;;    add the following lines to your .emacs file:
;;
;;       ; baci
;;       (autoload 'baci-mode "baci" "Enter baci mode." t)
;;       (setq auto-mode-alist (cons '("\\.dat\\'" . baci-mode)
;;                                   auto-mode-alist))
;;       (add-hook 'baci-mode-hook 'turn-on-font-lock)  ;; pbb
;;
;;    copy this file to Emacs site-lisp directory, or individual directory
;;    for instance ~/.site-lisp/baci.
;;
;;       cp baci.el [path to Emacs site-lisp directory] 
;;
;;    If you choose the latter, you need to provide in your .emacs file 
;;    the following line prior to the other lines above
;;
;;       (add-to-list 'load-path "~/.site-lisp/baci")
;;

(defvar baci-mode-hook nil)
(defvar baci-load-hook nil)

(defvar baci-font-lock-keywords
  (list
   '("^//.*" . font-lock-comment-face)  ;; comments
   '("[ ]+//.*" . font-lock-comment-face)  ;; comments
   '("^[-]+[-A-Z0-9 _&]+" . font-lock-reference-face)  ;; sections
   '("[/]+[A-Z0-9 _&'-]+" . font-lock-type-face)  ;; subsections
;   '("^[A-Z_]+" . font-lock-function-name-face) ;; core keywords
;   '("[A-Z_]+" . font-lock-keyword-face) ;; core keywords
;;   '("[\t]+" . highlight) ;; tabs
;;   '("[ ]+$" . highlight) ;; spaces just before end of line
  )
)

; imenu support stolen from python-mode :)

(defvar baci-imenu-generic-regexp "^--+\\([A-Za-z0-9 _&/]+\\)")

(defun baci-imenu-create-index-function ()
  "Function for finding Imenu definitions in baci.

Finds all definitions (classes, methods, or functions) in a baci
file for the Imenu package.

Returns a possibly nested alist of the form

        (INDEX-NAME . INDEX-POSITION)

The second element of the alist may be an alist, producing a nested
list as in

        (INDEX-NAME . INDEX-ALIST)
"
  (let (index-alist
        looking-p
        def-name
        def-pos
        )
    (goto-char (point-min))
    (setq looking-p
          (re-search-forward baci-imenu-generic-regexp (point-max) t))
    (while looking-p
      (setq def-pos (match-beginning 1))
      (setq def-name
            (buffer-substring-no-properties (match-beginning 1)
                                            (match-end 1)))
      (push (cons def-name def-pos) index-alist)
      (setq looking-p
            (re-search-forward baci-imenu-generic-regexp (point-max) t))
      )
    (nreverse index-alist)))

(defun baci-mode ()
  "Major mode for editing baci files."
  (interactive)
  (setq mode-name "baci")
  (setq major-mode 'baci-mode)
  ;; Install Imenu if available
  (when (baci-safe (require 'imenu))
    (setq imenu-create-index-function #'baci-imenu-create-index-function)
    ;;(setq imenu-generic-expression baci-imenu-generic-expression)
    (if (fboundp 'imenu-add-to-menubar)
        (imenu-add-to-menubar (format "%s-%s" "IM" mode-name)))
    )
  (make-local-variable 'font-lock-defaults)
  (setq font-lock-defaults '(baci-font-lock-keywords nil t))
  (run-hooks 'baci-mode-hook)
)

;; Utilities
(defmacro baci-safe (&rest body)
  "Safely execute BODY, return nil if an error occurred."
  `(condition-case nil
       (progn ,@ body)
     (error nil)))

(provide 'baci-mode)
(run-hooks 'baci-load-hook)


